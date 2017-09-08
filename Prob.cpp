#include <iomanip>
#include <Castro.H>
#include <Gravity.H>
#include "Castro_F.H"
#include <Gravity_F.H>

#define MAX_LEV 15

using namespace amrex;
namespace {

    // 2D->3D
    /*
    Real r_sync =  9.375000000000000e+07; // 3e9/32
    Real m_sync =  3.390985052132327e+33;
    Real g_sync = -2.575062324763961e+10;
    */

    /*
    Real r_sync =  4.687500000000000e+07; // 3e9/64
    Real m_sync =  3.386180724386528e+33;
    Real g_sync = -1.028565596622326e+11;
    */

    // 1D->3D
    Real r_sync =  1.875000000000000e+08; // 3e9/16
    Real m_sync =  2.690350324443699e+33;
    Real g_sync = -5.107527498930657e+09;

    /*
    Real r_sync =  9.375000000000000e+07; // 3e9/32
    Real m_sync =  2.689498567080453e+33;
    Real g_sync = -2.042364188022753e+10;
    */

    /*
    Real r_sync =  4.687500000000000e+07; // 3e9/64
    Real m_sync =  2.689274802310643e+33;
    Real g_sync = -8.168777057878873e+10;
    */

    /*
    Real r_sync =  2.343750000000000e+07; // 3e9/128
    Real m_sync =  2.689173263311251e+33;
    Real g_sync = -3.267387451683820e+11;
    */
}

void
Castro::problem_post_init ()
{
    // Note that restart does not call this function

    if (level > 0) return;

    const Geometry& geom = parent->Geom(level);
    const Real* dx   = geom.CellSize();
    int drdxfac = 1;//gravity->drdxfac;
//    int n1d = drdxfac * ( get_numpts() + 2 * NUM_GROW );
    int n1d = drdxfac * get_numpts();
    Real dr = dx[0] / double(drdxfac);

    const MultiFab& state = get_new_data(State_Type);
    
    Array<Real> radial_mass(n1d, 0.0);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    Array< Array<Real> > priv_radial_mass(nthreads);
    for (int i=0; i<nthreads; i++) {
		priv_radial_mass[i].resize(n1d,0.0);
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

    Array<Real> grav(n1d, 0.0);
    Array<Real> phi(n1d, 0.0);
    // Integrate radially outward to define the gravity
    ca_integrate_phi(radial_mass.dataPtr(),grav.dataPtr(),
                     phi.dataPtr(),&dr,&n1d);

    int index = r_sync/dr;
    BL_ASSERT(index < n1d);

    Real Gconst;
    get_grav_const(&Gconst);

    Real point_mass_g = -(g_sync - grav[index]) * r_sync * r_sync / Gconst;
    Real point_mass_m = m_sync - radial_mass[index];

    //set_pointmass(&point_mass);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << std::setprecision(15) << "problem_post_init:: r_sync       = " << r_sync << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: g_sync       = " << g_sync << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: g(r_sync)    = " << grav[index] << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: point_mass_g = " << point_mass_g << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: m_sync       = " << m_sync << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: m(r_sync)    = " << radial_mass[index] << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: point_mass_m = " << point_mass_m << std::endl;
    }

    
}

/*
void
Castro::problem_post_init ()
{
    int level = parent->finestLevel();
    int drdxfac = gravity->drdxfac;

    Array< Array<Real> > radial_mass(MAX_LEV);
    Array< Array<Real> > radial_vol(MAX_LEV);

    for (int lev = 0; lev <= level; lev++)
    {

//	const int NUM_STATE = getLevel(lev).get_new_data(State_Type).nComp();
	const int NUM_STATE = Castro::NUM_STATE;
        BoxArray ba = getLevel(lev).boxArray();

        // Create MultiFab with NUM_STATE components and no ghost cells
        MultiFab S(ba,NUM_STATE,0);
        S.copy(getLevel(lev).get_new_data(State_Type),0,0,NUM_STATE);

        if (lev < level)
        {
	    Castro* fine_level = dynamic_cast<Castro*>(&(parent->getLevel(lev+1)));
	    const MultiFab& mask = fine_level->build_fine_mask();
	    for (int n = 0; n < NUM_STATE; ++n)
		MultiFab::Multiply(S, mask, 0, n, 1, 0);
        }

        Castro& crse_level = dynamic_cast<Castro&>(parent->getLevel(lev));
        int n1d = drdxfac * crse_level.get_numpts();

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "problem_post_init:: lev = " << lev << std::endl;
            std::cout << "problem_post_init:: n1d = " << n1d << std::endl;
        }

        radial_mass[lev].resize(n1d);
        radial_vol[lev].resize(n1d);

        for (int i = 0; i < n1d; i++) radial_mass[lev][i] = 0.;
        for (int i = 0; i < n1d; i++) radial_vol[lev][i] = 0.;

        const Geometry& geom = parent->Geom(lev);
        const Real* dx   = geom.CellSize();
        Real dr = dx[0] / double(drdxfac);

#ifdef _OPENMP
	int nthreads = omp_get_max_threads();

	PArray< Array<Real> > priv_radial_mass(nthreads, PArrayManage);
	PArray< Array<Real> > priv_radial_vol (nthreads, PArrayManage);
	for (int i=0; i<nthreads; i++) {
	    priv_radial_mass.set(i, new Array<Real>(n1d,0.0));
	    priv_radial_vol.set (i, new Array<Real>(n1d,0.0));
	}
#pragma omp parallel
#endif
	{
#ifdef _OPENMP
	    int tid = omp_get_thread_num();
#endif
	    for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
	    {
	        const Box& bx = mfi.tilebox();
		FArrayBox& fab = S[mfi];

		ca_compute_radial_mass(bx.loVect(), bx.hiVect(), dx, &dr,
				       BL_TO_FORTRAN(fab),
#ifdef _OPENMP
				       priv_radial_mass[tid].dataPtr(),
				       priv_radial_vol[tid].dataPtr(),
#else
				       radial_mass[lev].dataPtr(),
				       radial_vol[lev].dataPtr(),
#endif
				       geom.ProbLo(),&n1d,&drdxfac,&lev);
	    }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
	    for (int i=0; i<n1d; i++) {
		for (int it=0; it<nthreads; it++) {
	            radial_mass[lev][i] += priv_radial_mass[it][i];
		    radial_vol [lev][i] += priv_radial_vol [it][i];
		}
	    }
#endif
	}

        ParallelDescriptor::ReduceRealSum(radial_mass[lev].dataPtr() ,n1d);
        ParallelDescriptor::ReduceRealSum(radial_vol[lev].dataPtr()  ,n1d);

    }

    int n1d = drdxfac * get_numpts();
    Array<Real> radial_mass_summed(n1d,0.0);

    // First add the contribution from this level
    for (int i = 0; i < n1d; i++)
    {
        radial_mass_summed[i] = radial_mass[level][i];
    }

    // Now add the contribution from coarser levels
    if (level > 0)
    {
        int ratio = parent->refRatio(level-1)[0];
        for (int lev = level-1; lev >= 0; lev--)
        {
            if (lev < level-1) ratio *= parent->refRatio(lev)[0];
            for (int i = 0; i < n1d/ratio; i++)
            {
                for (int n = 0; n < ratio; n++)
                {
                   radial_mass_summed[ratio*i+n] += 1./double(ratio) * radial_mass[lev][i];
                }
            }
        }
    }

    const Geometry& geom = parent->Geom(level);
    const Real* dx = geom.CellSize();
    Real dr        = dx[0] / double(drdxfac);

    // ***************************************************************** //
    // Compute the average density to use at the radius above
    //   max_radius_all_in_domain so we effectively count mass outside
    //   the domain.
    // ***************************************************************** //

    Array<Real> radial_vol_summed(n1d,0.0);
    Array<Real> radial_den_summed(n1d,0.0);

    // First add the contribution from this level
    for (int i = 0; i < n1d; i++)
         radial_vol_summed[i] =  radial_vol[level][i];

    // Now add the contribution from coarser levels
    if (level > 0)
    {
        int ratio = parent->refRatio(level-1)[0];
        for (int lev = level-1; lev >= 0; lev--)
        {
            if (lev < level-1) ratio *= parent->refRatio(lev)[0];
            for (int i = 0; i < n1d/ratio; i++)
            {
                for (int n = 0; n < ratio; n++)
                {
                   radial_vol_summed[ratio*i+n]  += 1./double(ratio) * radial_vol[lev][i];
                }
            }
        }
    }

    for (int i = 0; i < n1d; i++)
    {
        radial_den_summed[i] = radial_mass_summed[i];
        if (radial_vol_summed[i] > 0.) radial_den_summed[i]  /= radial_vol_summed[i];
    }

    // Compute the maximum radius at which all the mass at that radius is in the domain,
    //   assuming that the "hi" side of the domain is away from the center.
    Real max_radius_all_in_domain;
#if (BL_SPACEDIM > 1)
    Real center[3];
    get_center(center);
    Real x = Geometry::ProbHi(0) - center[0];
    Real y = Geometry::ProbHi(1) - center[1];
    max_radius_all_in_domain = std::min(x,y);
#if (BL_SPACEDIM == 3)
    Real z = Geometry::ProbHi(2) - center[2];
    max_radius_all_in_domain = std::min(max_radius_all_in_domain,z);
#endif
#endif

    Array<Real> radial_grav(n1d, 0.0);


    // Integrate radially outward to define the gravity
    ca_integrate_grav(radial_mass_summed.dataPtr(),radial_den_summed.dataPtr(),
		      radial_grav.dataPtr(),&max_radius_all_in_domain,&dr,&n1d);

    int index = r_sync/dr;
    BL_ASSERT(index < n1d);

    Real Gconst;
    get_grav_const(&Gconst);

    Real point_mass_g = -(g_sync - radial_grav[index]) * r_sync * r_sync / Gconst;
    Real point_mass_m = m_sync - radial_mass_summed[index];

    //set_pointmass(&point_mass);

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << std::setprecision(15) << "problem_post_init:: r_sync       = " << r_sync << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: g_sync       = " << g_sync << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: g(r_sync)    = " << radial_grav[index] << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: point_mass_g = " << point_mass_g << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: m_sync       = " << m_sync << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: m(r_sync)    = " << radial_mass_summed[index] << std::endl;
	std::cout << std::setprecision(15) << "problem_post_init:: point_mass_m = " << point_mass_m << std::endl;
    }

}
*/
