#include <cstdlib>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cux/cux.h"


static int cuda_initialized = 0;	// CUDA initialized
static int cuda_enabled = 0;		// Should GPU acceleration be used?

int main(int argc, char **argv)
{
	int dev;
	const char *devStr = getenv("CUDA_DEVICE");
	bool autoselect = devStr == NULL;

	if(!autoselect)
	{
		dev = atoi(devStr);

		// disable GPU acceleration
		if(dev == -1)
		{
			cuda_initialized = 1;
			cuda_enabled = 0;

			std::cerr << "GPU accelerator: Using CPU\n";
			return true;
		}
		else
		{
			cuxErrCheck( cudaSetDevice(dev) );
		}
	}


        // get device properties
        cudaDeviceProp deviceProp;
        cuxErrCheck( cudaGetDeviceProperties(&deviceProp, dev) );

	char buf[1000];
	snprintf(buf, sizeof(buf), "GPU accelerator: Using Device %d: \"%s\"%s\n", dev, deviceProp.name, autoselect ? " (autoselected)" : "");
	std::cerr << buf;
};


