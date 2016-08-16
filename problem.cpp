#include <Castro.H>
#include "Castro_F.H"
#include <Gravity_F.H>

namespace {
    Real r_sync = 0.5e9;
    Real m_sync = 1.e33;  // Enclosed mass inside r_sync
}

void
Castro::problem_post_init ()
{
    // Note that restart does not call this function

    if (level > 0) return;

    int n1d = get_numpts();
    int drdxfac = 1;
    const Geometry& geom = parent->Geom(level);
    const Real* dx   = geom.CellSize();
    Real dr = dx[0] / double(drdxfac);

    const MultiFab& state = get_new_data(State_Type);
    
    Array<Real> radial_mass(n1d, 0.0);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    PArray< Array<Real> > priv_radial_mass(nthreads, PArrayManage);
    for (int i=0; i<nthreads; i++) {
	priv_radial_mass.set(i, new Array<Real>(n1d,0.0));
    }
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
	int tid = omp_get_thread_num();
	Array<Real>& mass = priv_radial_mass[tid];
#else
	Array<Real>& mass = radial_mass;
#endif
	
	Array<Real> vol(n1d, 0.0);
	
	for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();
	    ca_compute_radial_mass(bx.loVect(), bx.hiVect(), dx, &dr,
				   BL_TO_FORTRAN_N(state[mfi],Density),
				   mass.dataPtr(),
				   vol.dataPtr(),
				   geom.ProbLo(),&n1d,&drdxfac,&level);
	}

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
	for (int i=0; i<n1d; i++) {
	    for (int it=0; it<nthreads; it++) {
		radial_mass[i] += priv_radial_mass[it][i];
	    }
	}
#endif
    }

    ParallelDescriptor::ReduceRealSum(radial_mass.dataPtr(),n1d);

    int index = r_sync/dr;
    BL_ASSERT(index < n1d);

    Castro::point_mass = m_sync - radial_mass[index];

    set_pointmass(&point_mass);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "problem_post_init:: reset pointmass to " << point_mass
		  << ".  Do NOT forget to update inputs file." << std::endl;
    }

    
}

