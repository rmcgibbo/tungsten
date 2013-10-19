#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <streambuf>
#include <sstream>

#include "mpi.h"
#include "OpenMM.h"
#include "openmm/serialization/XmlSerializer.h"

#define MASTER 0
#define PLATFORM_NAME "OpenCL"

// Prototypes
int main(int argc, char* argv[]);
OpenMM::Context* createContext(const char* systemXml, const char* integratorXml);
void setStateFromXML(OpenMM::Context& context, const char* stateXML);


int main(int argc, char* argv[]) {

  int rank;       // id of this node
  int n_proc;     // number of processors
  int error = 0;    // has an error occured?
  int system_xml_length = 0;
  char* system_xml;
  int integrator_xml_length = 0;
  char* integrator_xml;
  int state_xml_length = 0;
  char* state_xml;
  MPI::Init(argc, argv);

  n_proc = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
  
  // Parse the command line
  if (rank == MASTER) {
    if (argc != 4) {
      fprintf(stderr, "usage: %s <system.xml> <integrator.xml> <state.xml>\n", argv[0]);
      error = 1;
    } else { 
      // Read system_xml file
      FILE* fp = fopen(argv[1], "r");
      fseek(fp, 0, SEEK_END);
      system_xml_length = ftell(fp);
      fseek(fp, 0, SEEK_SET);
      system_xml = (char*) malloc(system_xml_length*sizeof(char));
      fread(system_xml, 1, system_xml_length, fp);
      fclose(fp);

      fp = fopen(argv[2], "r");
      fseek(fp, 0, SEEK_END);
      integrator_xml_length = ftell(fp);
      fseek(fp, 0, SEEK_SET);
      integrator_xml = (char*) malloc(integrator_xml_length*sizeof(char));
      fread(integrator_xml, 1, integrator_xml_length, fp);
      fclose(fp);

      fp = fopen(argv[3], "r");
      fseek(fp, 0, SEEK_END);
      state_xml_length = ftell(fp);
      fseek(fp, 0, SEEK_SET);
      state_xml = (char*) malloc(state_xml_length*sizeof(char));
      fread(state_xml, 1, state_xml_length, fp);
      fclose(fp);
    }
  }
  // Check for any configuration errors from command line parsing.
  MPI::COMM_WORLD.Bcast(&error, 1, MPI_INT, MASTER);
  if (error != 0) {
    if (rank == MASTER) {
      fprintf(stderr, "Error: Program terminated with error code %d\n", error);
    }
    MPI::Finalize();
    exit(error);
  }
 
  // Broadcast the contents of the system and integator xml files
  MPI::COMM_WORLD.Bcast(&system_xml_length, 1, MPI_INT, MASTER);
  MPI::COMM_WORLD.Bcast(&integrator_xml_length, 1, MPI_INT, MASTER);
  MPI::COMM_WORLD.Bcast(&state_xml_length, 1, MPI_INT, MASTER);
  if (rank != MASTER) {
    system_xml = (char*) malloc(system_xml_length);
    integrator_xml = (char*) malloc(integrator_xml_length);
    state_xml = (char*) malloc(state_xml_length);
  }
  MPI::COMM_WORLD.Bcast(system_xml, system_xml_length + 1, MPI_CHAR, MASTER);
  MPI::COMM_WORLD.Bcast(integrator_xml, integrator_xml_length + 1, MPI_CHAR, MASTER);
  MPI::COMM_WORLD.Bcast(state_xml, state_xml_length + 1, MPI_CHAR, MASTER);

  // Build the context and set the state
  OpenMM::Context* context =  createContext(system_xml, integrator_xml);
  setStateFromXML(*context, state_xml);
  free(system_xml);
  free(integrator_xml);
  free(state_xml);


  OpenMM::Integrator* integrator = &(context->getIntegrator());
  // MAIN LOOP
  for (int i = 0; i < 1; i++)
    integrator->step(100);



  
  printf("I am rank=%d\n", rank);


  delete context;
  MPI::Finalize();
}


OpenMM::Context* createContext(const char* systemXml, const char* integratorXml) {
  int rank = MPI::COMM_WORLD.Get_rank();
  std::stringstream systemXmlStream;
  std::stringstream integratorXmlStream;
  systemXmlStream << systemXml;
  integratorXmlStream << integratorXml;

  OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

  OpenMM::System* system = OpenMM::XmlSerializer::deserialize<OpenMM::System>(systemXmlStream);
  OpenMM::Integrator* integrator = OpenMM::XmlSerializer::deserialize<OpenMM::Integrator>(integratorXmlStream);
  OpenMM::Platform* platform = &OpenMM::Platform::getPlatformByName(PLATFORM_NAME);


  OpenMM::Context* context = new OpenMM::Context(*system, *integrator, *platform);
  if (rank == MASTER) {
    printf("Context created on OpenMM platform: %s\n", (*context).getPlatform().getName().c_str());
    std::vector<std::string> names = (*context).getPlatform().getPropertyNames();
    for (int i = 0; i < names.size(); i++) {
      std::string value = (*context).getPlatform().getPropertyValue(*context, names[i]);
      printf("  %s: %s\n", names[i].c_str(), value.c_str());
    }
  }

  return context;
}


void setStateFromXML(OpenMM::Context& context, const char* stateXML) {
  std::stringstream s;
  s << stateXML;
  OpenMM::State* state = OpenMM::XmlSerializer::deserialize<OpenMM::State>(s);
  context.setState(*state);
}
