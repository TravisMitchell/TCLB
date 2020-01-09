#ifndef CBSTOPMAX_H
#define CBSTOPMAX_H

#include "../CommonHandler.h"

#include "vHandler.h"
#include "Callback.h"

class  cbStopMax  : public  Callback  {
        std::vector< int > what;
        std::vector< double > MAX;
	int times, score;
	int old_iter_type;
	public:
	static std::string xmlname;
int Init ();
int DoIt ();
int Finish ();
};

#endif // CBSTOPMAX_H
