/*! \file hermite_adap_gpu_integrator_body.cu
 *  \brief inluded multiple times to acheive template specialization for various precisions
 *
 *  This is included from hermite_adap_gpu.cu like
 *  template<>
 *  void gpu_hermite_adap_integrator<double,double>::integrate(gpu_ensemble &ens, double dT)
 *  #include"hermite_adap_gpu_integrator_body.cu"
 *
*/

//template<>
//void gpu_hermite_adap_integrator<real_hi,real_lo>::integrate(gpu_ensemble &ens, double dT)
{ 
	/* Upload the kernel parameters */ 
	if(ens.last_integrator() != this) 
	{ 
		ens.set_last_integrator(this); 
		configure_grid(gridDim, threadsPerBlock, ens.nsys()); 
		cudaMemcpyToSymbol(gpu_hermite_adap_ens,
				   &ens, 
				   sizeof(gpu_hermite_adap_ens) ); 
		if(dT == 0.) { return; } 
	} 
 
	// flush CPU/GPU output logs
	log::flush(log::memory | log::if_full);

	if(ens.nbod()==3) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,3><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,3><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,3><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else if(ens.nbod()==4) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,4><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,4><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,4><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else if(ens.nbod()==5) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,5><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,5><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,5><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else if(ens.nbod()==6) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,6><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,6><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,6><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else if(ens.nbod()==7) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,7><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,7><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,7><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else if(ens.nbod()==8) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,8><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,8><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,8><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else if(ens.nbod()==9) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,9><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,9><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,9><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else if(ens.nbod()==10) {
		// execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,10><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,10><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,10><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else {
#ifdef __LARGE_N_SYSTEM__
        if (ens.nbod()==N_BODY)  {
                // execute the kernel
		switch(prec){
			// double precision
			case 1:
				gpu_hermite_adap_integrator_kernel<1,N_BODY><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// signle precision
			case 2:
				gpu_hermite_adap_integrator_kernel<2,N_BODY><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
				// mixed precision
			case 3:
				gpu_hermite_adap_integrator_kernel<3,N_BODY><<<gridDim, threadsPerBlock>>>(dT, h, stepfac);
				break;
		}
	}
	else {
		// How do we get an error message out of here?
		ERROR("Invalid number of bodies using definition N_BODY. Ensemble should equal N_BODY");
		return;
	}
#endif
	// How do we get an error message out of here?
	ERROR("Invalid number of bodies. (only up to 10 bodies per system)");
	return;
	}
	//printf("%s\n", cudaGetErrorString(cudaGetLastError())); 

	// flush CPU/GPU output logs
	log::flush(log::memory);
};

