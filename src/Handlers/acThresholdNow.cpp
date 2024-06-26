#include "acThresholdNow.h"

std::string acThresholdNow::xmlname = "ThresholdNow";

int acThresholdNow::Init () {
		double val;
		double * start = NULL;
		double * slice = NULL;
		GenericAction::Init();
		pugi::xml_attribute attr = node.attribute("Level");
		if (attr) {
			level = attr.as_int();
		} else {
		        level = 0.5;
		}
		DEBUG_M;
//		par = solver->getPars();
		par = NumberOfParameters();
		DEBUG_M;
		if (solver->mpi_rank == 0) {
        		if (par == 0) {
        			ERROR("No parameters defined!\n");
        			return -1;
        		}			
			output("Parameters: %d\n", par);
			start = new double[par];
			slice = new double[par];
		}
		DEBUG_M;
//		solver->getPar(start);
		GetParameters(start);
		DEBUG_M;
		int msg=0;
		const Model::Setting& it = solver->lattice->model->settings.by_name("Threshold");
		if (!it) {
			ERROR("'Threshold' is not a setting");
			return -1;
		}
	        solver->lattice->SetSetting(it, level);
	        if (slice != NULL) for (int j=0;j<par;j++) slice[j]=start[j]>level ? 1.0 : 0.0;
//       		solver->setPar(slice);
		SetParameters(slice);
		return 0;
	}


// Register the handler (basing on xmlname) in the Handler Factory
template class HandlerFactory::Register< GenericAsk< acThresholdNow > >;
