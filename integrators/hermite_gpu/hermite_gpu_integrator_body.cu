//template<>
//void gpu_hermite_integrator<real_hi,real_lo>::integrate(gpu_ensemble &ens, double dT)
{ 
	/* Upload the kernel parameters */ 
	if(ens.last_integrator() != this) 
	{ 
		ens.set_last_integrator(this); 
		configure_grid(gridDim, threadsPerBlock, ens.nsys()); 
		cudaMemcpyToSymbol(gpu_hermite_ens,
				   &ens, 
				   sizeof(gpu_hermite_ens) ); 
		if(dT == 0.) { return; } 
	} 
 
	if(ens.nbod()==3) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_integrator_kernel<1,3><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
				// signle precision
			case 2:
				gpu_hermite_integrator_kernel<2,3><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
				// mixed precision
			case 3:
				gpu_hermite_integrator_kernel<3,3><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
		}
	}
	else if(ens.nbod()==4) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_integrator_kernel<1,4><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
				// signle precision
			case 2:
				gpu_hermite_integrator_kernel<2,4><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
				// mixed precision
			case 3:
				gpu_hermite_integrator_kernel<3,4><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
		}
	}
	else if(ens.nbod()==5) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_integrator_kernel<1,5><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
				// signle precision
			case 2:
				gpu_hermite_integrator_kernel<2,5><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
				// mixed precision
			case 3:
				gpu_hermite_integrator_kernel<3,5><<<gridDim, threadsPerBlock>>>(dT, h);
				break;
		}
	}
	else {
	// How do we get an error message out of here?
	ERROR("Invalid number of bodies.");
	return;
	}



	printf("%s\n", cudaGetErrorString(cudaGetLastError())); 
};

