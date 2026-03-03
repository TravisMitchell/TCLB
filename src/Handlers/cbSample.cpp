#include "cbSample.h"
std::string cbSample::xmlname = "Sample";
#include "../HandlerFactory.h"

int cbSample::Init () {
		std::string nm="Sampler";
		Callback::Init();
		if (everyIter == 0) {
                	error("Iteration value in sampler should not be zero");
                	return -1;
                }
		pugi::xml_attribute attr=node.attribute("what");
		if (attr) {
			s.add_from_string(attr.value(),',');
		}
		else {
			s.add_from_string("all",',');
		}
		for (pugi::xml_node par = node.first_child(); par; par = par.next_sibling()) {
			if (strcmp(par.name(),"Point") == 0) {
				lbRegion loc;
				attr = par.attribute("dx");
				if (attr) {
					loc.dx = solver->units.alt(attr.value());
				}
				attr = par.attribute("dy");
				if (attr) {
					loc.dy = solver->units.alt(attr.value());
				}
				attr = par.attribute("dz");
				if (attr) {
					loc.dz = solver->units.alt(attr.value());
				}
				loc = solver->lattice->getLocalBoundingBox().intersect(loc);

				if (loc.nx == 1) {
					unsigned int lid = 0;
					auto variant = solver->getLatticeVariant();
					if (auto* lattice = std::get_if<Lattice<ArbLattice>*>(&variant)) {
						// cache lid for arbitrary lattice
						const real_t offset = 0.5;
						vector_t point{real_t(loc.dx) + offset, real_t(loc.dy) + offset, real_t(loc.dz) + offset};
						lid = (*lattice)->getCartesianCoordinateLid(point);
					}
					solver->lattice->sample->addPoint(loc, solver->mpi_rank, lid);
				}

			} else {
				error("Uknown element in Sampler\n");
				return -1;
			}
		}
		filename = solver->outIterFile(nm, ".csv");
		solver->lattice->sample->units = &solver->units;
		solver->lattice->sample->mpi_rank = solver->mpi_rank;
		solver->lattice->sample->Allocate(&s,startIter,everyIter);
		solver->lattice->sample->initCSV(filename.c_str());
		return 0;
		}


int cbSample::DoIt () {
		Callback::DoIt();
		solver->lattice->sample->writeHistory(solver->iter);
		solver->lattice->sample->startIter = solver->iter;
		return 0;
		}


int cbSample::Finish () {
	   solver->lattice->sample->Finish();
	   return Callback::Finish();
	 }	 


// Register the handler (basing on xmlname) in the Handler Factory
template class HandlerFactory::Register< GenericAsk< cbSample > >;
