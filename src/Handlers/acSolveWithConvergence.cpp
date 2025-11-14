#include "acSolveWithConvergence.h"
#include <climits>

std::string acSolveWithConvergence::xmlname = "SolveWithConvergence";

int acSolveWithConvergence::Init () {
	GenericAction::Init();
	pugi::xml_attribute attr;
	
	// Get the main action name (run once per outer iteration)
	attr = node.attribute("action");
	std::string action_name;
	int action_id = -1;
	if (attr) {
		action_name = attr.value();
		const Model::Action& act = solver->lattice->model->actions.by_name(action_name);
		if (!act) {
			ERROR("Unknown Action: %s", action_name.c_str());
			return -1; 	
		}
		action_id = act.id;
	}

	// Get the convergence action name (looped until convergence)
	attr = node.attribute("converge_action");
	std::string converge_action_name;
	int converge_action_id = -1;
	if (attr) {
		converge_action_name = attr.value();
		const Model::Action& act = solver->lattice->model->actions.by_name(converge_action_name);
		if (!act) {
			ERROR("Unknown Convergence Action: %s", converge_action_name.c_str());
			return -1; 	
		}
		converge_action_id = act.id;
	} else {
		ERROR("Have to specify the convergence action name in SolveWithConvergence (attribute 'converge_action')");
		return -1;
	}		// Get the convergence global name
		attr = node.attribute("global");
		std::string global_name;
		if (attr) {
			global_name = attr.value();
		}
		if (global_name == "") {
			ERROR("Have to specify the global name for convergence checking in SolveWithConvergence");
			return -1;
		}
		int global_index = solver->getGlobal(global_name.c_str());
		if (global_index < 0) {
			ERROR("Unknown Global: %s", global_name.c_str());
			return -1;
		}
		
		// Get the convergence threshold
		attr = node.attribute("threshold");
		double threshold = 1e-6;  // default
		if (attr) {
			threshold = solver->units.alt(attr.value());
		}
		
		// Get the check frequency (every N iterations)
		attr = node.attribute("check_every");
		int check_every = 1;  // default: check every iteration
		if (attr) {
			check_every = atoi(attr.value());
			if (check_every < 1) check_every = 1;
		}
		
		// Get minimum iterations before checking convergence
		attr = node.attribute("min_iterations");
		int min_iterations = 10;  // default
		if (attr) {
			min_iterations = atoi(attr.value());
		}
		
		// Get maximum iterations (safety limit)
		attr = node.attribute("max_stage_iterations");
		int max_stage_iterations = INT_MAX;  // default: no limit on stage iterations
		if (attr) {
			max_stage_iterations = atoi(attr.value());
			if (max_stage_iterations < 1) max_stage_iterations = INT_MAX;
		}

		// Get outer iterations count
		attr = node.attribute("iterations");
		int iterations = 1; // default: perform one outer iteration
		if (attr) {
			iterations = atoi(attr.value());
			if (iterations < 1) iterations = 1;
		}
		
	
	if (solver->mpi_rank == 0) {
		output("SolveWithConvergence: action=%s, converge_action=%s, global=%s, threshold=%g, check_every=%d, min_iter=%d, iterations=%d, max_stage_iter=%d\n",
		       (action_id >= 0 ? action_name.c_str() : "<none>"),
		       converge_action_name.c_str(),
		       global_name.c_str(), threshold, check_every, min_iterations, iterations, max_stage_iterations);
	}		// Execute internal initialization (sets up timing, etc.)
		if (GenericAction::ExecuteInternal()) return -1;
		
		// Allocate array for globals
		int n_globals = solver->lattice->model->globals.size();
		real_t * glob = new real_t[n_globals];
		
		int stop = 0;
		int converged = 0;
		
	for (int outer = 0; outer < iterations && !stop; ++outer) {
		// Advance the global iteration counter once per outer iteration
		solver->steps = 1;
		solver->iter += solver->steps;
		
		// Run the main action once if specified
		if (action_id >= 0) {
			if (solver->mpi_rank == 0) {
				output("  Running main action: %s\n", action_name.c_str());
			}
			solver->lattice->IterateAction(action_id, 1, solver->iter_type);
			CudaDeviceSynchronize();
			MPI_Barrier(MPMD.local);
		}
		
		// Reset convergence iteration counter for this outer iteration
		int stage_iter_count = 0;
		converged = 0;
		
		// Inner loop: iterate the convergence action until converged or max_stage_iterations
		while (!converged && stage_iter_count < max_stage_iterations) {
			// Number of action iterations to perform in this batch
			int stage_steps = 1;  // Do one iteration at a time for precise convergence checking
			
			stage_iter_count += stage_steps;
			
			// Iterate the convergence action with global calculation
			// Set ITER_GLOBS to calculate globals after iteration
			// NOTE: We don't increment solver->iter here - only the outer loop does
			int iter_type_with_globs = solver->iter_type | ITER_GLOBS | ITER_LASTGLOB;
			solver->lattice->IterateAction(converge_action_id, stage_steps, iter_type_with_globs);
			
			CudaDeviceSynchronize();
			MPI_Barrier(MPMD.local);
			
			// Check convergence every check_every iterations
			if ((stage_iter_count >= min_iterations) && (stage_iter_count % check_every == 0)) {
				// Get globals
				solver->lattice->getGlobals(glob);
				
				// Get the specific convergence metric
				double conv_value = glob[global_index];
				
				// Broadcast convergence value to all ranks
				MPI_Bcast(&conv_value, 1, MPI_DOUBLE, 0, MPMD.local);
				
				if (solver->mpi_rank == 0) {
					output("  Convergence action iter %d (total %d): %s = %g (threshold = %g)\n", 
					       stage_iter_count, solver->iter, global_name.c_str(), conv_value, threshold);
				}
				
				// Check if converged
				if (fabs(conv_value) < threshold) {
					converged = 1;
					if (solver->mpi_rank == 0) {
						notice("Convergence action converged after %d iterations: %s = %g < %g\n",
						       stage_iter_count, global_name.c_str(), conv_value, threshold);
					}
				}
			}
		}
		
		// Warn if convergence action didn't converge within max_stage_iterations
		if (!converged && stage_iter_count >= max_stage_iterations) {
			if (solver->mpi_rank == 0) {
				warning("Convergence action reached max_stage_iterations (%d) without convergence, continuing to next outer iteration\n", 
				        max_stage_iterations);
			}
		}
		
		// Execute callbacks at outer iteration boundaries
		for (size_t i=0; i<solver->hands.size(); i++) {
			if (solver->hands[i].Now(solver->iter)) {
				int ret = solver->hands[i].DoIt();
				switch (ret) {
				case ITERATION_STOP:
					stop=1;
				case 0:
					break;
				default:
					delete[] glob;
					return -1;
				}
			}
		}
		
		if (stop) break;
	}
		
		delete[] glob;
		
		CudaDeviceSynchronize();
		MPI_Barrier(MPMD.local);
		GenericAction::Unstack();
		return 0;
	}


// Register the handler (basing on xmlname) in the Handler Factory
template class HandlerFactory::Register< GenericAsk< acSolveWithConvergence > >;
