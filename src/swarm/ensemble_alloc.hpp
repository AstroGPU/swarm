#pragma once
#include "types/ensemble.hpp"

namespace swarm {

/*! Allocator based version of ensemble containing memory management routines
 * It takes an allocator as a parameter and uses the allocator for 
 * allocate, deallocate and copying the ensemble. The allocated memories are
 * protected by shared pointers.
 */
template< int W , template<class T> class _Allocator >
struct EnsembleAlloc : public EnsembleBase<W> {
        typedef EnsembleBase<W> Base;
        typedef EnsembleAlloc Self;
        typedef typename Base::Body Body;
        typedef typename Base::Sys Sys;
        typedef _Allocator<Body> BodyAllocator;
        typedef _Allocator<Sys> SysAllocator;

        typedef boost::shared_ptr<Body> PBody;
        typedef boost::shared_ptr<Sys> PSys;

        //! Deep copy of this ensemble
        EnsembleAlloc clone() {
                return cloneTo<EnsembleAlloc>();
        }

        //! Create a new ensemble that can accomodate nsys systems with nbod bodies
        //! Arrays are allocated on the heap but ensemble structure is value-copied
        static EnsembleAlloc create(const int& nbod, const int& nsys) {
                PBody b ( BodyAllocator::alloc( Base::body_element_count(nbod,nsys) ), &BodyAllocator::free );
                PSys s ( SysAllocator::alloc( Base::sys_element_count(nsys) ), &SysAllocator::free );
                return EnsembleAlloc(nbod,nsys,b,s);
        }

        //! Clone to a different memory (e.g. GPU)
        template< class Other > 
        Other cloneTo() {
                Other o = Other::create(Base::nbod(),Base::nsys());
                copyTo(o);
                return o;
        }

        //! Copy to another ensemble located on a different memory
        template< class Other > 
        void copyTo(Other& o){
                assert(o.nsys() == Base::nsys() && o.nbod() == Base::nbod());
                alloc_copy(BodyAllocator(),typename Other::BodyAllocator(),Base::bodies().begin(),Base::bodies().end(),o.bodies().begin());
                alloc_copy(SysAllocator(),typename Other::SysAllocator(),Base::systems().begin(),Base::systems().end(),o.systems().begin());
        }

        EnsembleAlloc(){}
        private:
        EnsembleAlloc(const int& nbod,const int& nsys, PBody b, PSys s):
                Base(nbod,nsys,b.get(),s.get())
                ,_body(b),_sys(s){}

        PSys _sys;
        PBody _body;
};


//! Default ensemble class for most of uses
typedef EnsembleAlloc< ENSEMBLE_CHUNK_SIZE , DefaultAllocator > defaultEnsemble;
//! Ensemble allocated on host memory
typedef EnsembleAlloc< ENSEMBLE_CHUNK_SIZE , DefaultAllocator > hostEnsemble;
//! Ensemble allocated on [GPU] device memory
typedef EnsembleAlloc< ENSEMBLE_CHUNK_SIZE , DeviceAllocator > deviceEnsemble;


//! Provided for backward compatibility
typedef hostEnsemble cpu_ensemble;
//! Provided for backward compatibility
typedef deviceEnsemble gpu_ensemble;

}