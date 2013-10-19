a.out: mainloop.cpp
	mpic++ mainloop.cpp -I$(OPENMM_INCLUDE_PATH) -L$(OPENMM_LIB_PATH) -lOpenMM -lOpenMMSerialization

clean:
	rm a.out
