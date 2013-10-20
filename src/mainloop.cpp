#include <cstdlib>
#include <sys/stat.h>
#include <cstdio>
#include <string>
#include <fstream>
#include <sstream>

#include "mpi.h"
#include "INIReader.h"
#include "OpenMM.h"
#include "openmm/serialization/XmlSerializer.h"

#define MASTER 0
//#define PLATFORM_NAME "OpenCL"
#define PLATFORM_NAME "Reference"

// ************************************************************************* //
// Prototypes
// ************************************************************************* //
typedef struct{
    int n_rounds;
    int n_steps_per_round;
    int save_frequency;
    std::string output_root_path;
} ConfigOpts;


/**
 * Main entry point
 */
int main(int argc, char* argv[]);

/**
 * Boot up a context from a serialized system and integrator
 */
OpenMM::Context* createContext(std::ifstream& systemXml, std::ifstream& integratorXml);

/**
 * Convenience method to exit from MPI
 */
void exitWithMessage(const char* message);

/**
 * Parse the config file
 */ 
void parseConfigFile(const char* configFileName, ConfigOpts* out);



// ************************************************************************* //
// Functions
// ************************************************************************* //
int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    int n_proc = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();

    // Parse the command line
    if (argc != 5) {
        fprintf(stderr, "usage: %s <system.xml> <integrator.xml> <state.xml> <config.ini>\n", argv[0]);
        if (rank == MASTER)
            fprintf(stderr, "n_rounds must be given and greater than 0\n");
        MPI::Finalize();
        exit(EXIT_FAILURE);
    }

    // Create the context from the input files
    std::ifstream systemXml(argv[1]);
    std::ifstream integratorXml(argv[2]);
    OpenMM::Context* context =  createContext(systemXml, integratorXml);
    OpenMM::Integrator& integrator = context->getIntegrator();

    // Set the initial state
    std::ifstream stateXml(argv[3]);
    OpenMM::State* state = OpenMM::XmlSerializer::deserialize<OpenMM::State>(stateXml);
    context->setState(*state);

    // Get the config file
    ConfigOpts opts;
    parseConfigFile(argv[4], &opts);
    
    // Create the output directory for our work on this node
    std::stringstream s;
    s <<  opts.output_root_path << "/node-" << rank;
    std::string outdir = s.str();
    mkdir(outdir.c_str(), 0777);
    
    for (int round = 0; round < opts.n_rounds; round++) {
        for (int step = 0; step < opts.n_steps_per_round; step += opts.save_frequency) {
            integrator.step(opts.save_frequency);
    
        }
    }
    

    delete context;
    MPI::Finalize();
}
//   // Build the context and set the state
//   OpenMM::Context* context =  createContext(system_xml, integrator_xml);
//   setStateFromXML(*context, state_xml);
//   free(system_xml);
//   free(integrator_xml);
//   free(state_xml);
//   OpenMM::Integrator* integrator = &(context->getIntegrator());
// 
// 
// 
//   
//   // MAIN LOOP
//   for (int i = 0; i < 1; i++)
//     integrator->step(10);
// 
// 
// 
//   
//   printf("I am rank=%d\n", rank);
// 
// 
//   delete context;
//   MPI::Finalize();
// }


OpenMM::Context* createContext(std::ifstream& systemXml, std::ifstream& integratorXml) {
  int rank = MPI::COMM_WORLD.Get_rank();
  OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

  OpenMM::System* system = OpenMM::XmlSerializer::deserialize<OpenMM::System>(systemXml);
  OpenMM::Integrator* integrator = OpenMM::XmlSerializer::deserialize<OpenMM::Integrator>(integratorXml);
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


void exitWithMessage(const char* message) {
  int rank = MPI::COMM_WORLD.Get_rank();
  if (rank == MASTER)
      fprintf(stderr, "%s\n", message);
  MPI::Finalize();
  exit(EXIT_FAILURE);
}


void parseConfigFile(const char* configFileName, ConfigOpts* out) {
  INIReader reader(configFileName);
  if (reader.ParseError() < 0) {
      exitWithMessage("Could not find config file");
  }
  
  out->n_steps_per_round = reader.GetInteger("", "n_steps_per_round", -1);
  out->n_rounds = reader.GetInteger("", "n_rounds", -1);
  out->save_frequency = reader.GetInteger("", "save_frequency", -1);
  out->output_root_path = reader.Get("", "output_root_path", ".");

  
  if (out->n_steps_per_round <= 0)
      exitWithMessage("n_steps_per_round must be given and greater than 0");
  if (out->n_rounds <= 0)
    exitWithMessage("n_rounds must be given and greater than 0");
  if (out->save_frequency <= 0)
    exitWithMessage("save_frequency must be given and greater than 0");
}