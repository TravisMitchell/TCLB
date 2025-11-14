#ifndef ACSOLVEWITHCONVERGENCE_H
#define ACSOLVEWITHCONVERGENCE_H

#include "../CommonHandler.h"

#include "vHandler.h"
#include "Action.h"
#include "GenericAction.h"

class  acSolveWithConvergence  : public  GenericAction  {
	public:
	static std::string xmlname;
int Init ();
};

#endif // ACSOLVEWITHCONVERGENCE_H
