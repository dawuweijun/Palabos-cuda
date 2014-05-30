/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* This code was written with help of a Fortran code kindly provided
 * by Prof. Taehun Lee, and is massively co-authored by
 * Andrea Parmigiani.
 */

#ifndef HE_LEE_PROCESSOR_2D_HH
#define HE_LEE_PROCESSOR_2D_HH

#include "heLeeProcessor2D.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
namespace plb
{

template<typename T>
class FiniteDifference2D
{
public:
    FiniteDifference2D ( ScalarField2D<T> const& field_,
                         plint x_, plint y_ )
        : field ( field_ ),
          x ( x_ ), y ( y_ )
    { }
    T val ( plint iX, plint iY ) const
    {
        return field.get ( iX,iY );
    }
    void centralGradient ( Array<T,2>& gradient ) const
    {
        const T c1 = ( T ) 1./ ( T ) 3.;
        const T c2 = ( T ) 1./ ( T ) 12.;
        gradient[0] =
            ( val ( x+1,y )   - val ( x-1,y ) ) *c1
            + ( val ( x+1,y+1 ) - val ( x-1,y-1 ) + val ( x+1,y-1 ) - val ( x-1,y+1 ) ) *c2;
        gradient[1] =
            ( val ( x,y+1 )   - val ( x,y-1 ) ) *c1
            + ( val ( x+1,y+1 ) - val ( x-1,y-1 ) + val ( x-1,y+1 ) - val ( x+1,y-1 ) ) *c2;
    }
    void biasedGradient ( Array<T,2>& gradient ) const
    {
        const T c1 = ( T ) 1./ ( T ) 6.;
        const T c2 = ( T ) 1./ ( T ) 24.;
        gradient[0] =
            ( -   val ( x+2,y )   + 4.*val ( x+1,y )
              - 4.*val ( x-1,y ) +   val ( x-2,y ) ) *c1+
            ( - val ( x+2,y+2 ) + 4.*val ( x+1,y+1 )
              -4.*val ( x-1,y-1 ) + val ( x-2,y-2 )
              - val ( x+2,y-2 ) + 4.*val ( x+1,y-1 )
              -4.*val ( x-1,y+1 ) + val ( x-2,y+2 ) ) *c2;
        gradient[1] =
            ( -  val ( x,y+2 ) + 4.*val ( x,y+1 )
              - 4.*val ( x,y-1 ) + val ( x,y-2 ) ) *c1
            + ( -val ( x+2,y+2 ) + 4.*val ( x+1,y+1 )
                -4.*val ( x-1,y-1 ) + val ( x-2,y-2 )
                -val ( x-2,y+2 ) + 4.*val ( x-1,y+1 )
                -4.*val ( x+1,y-1 ) + val ( x+2,y-2 ) ) *c2;
    }

    T laplacian() const
    {
        static const T c1 = ( T ) 2./ ( T ) 3.;
        static const T c2 = ( T ) 1./ ( T ) 6.;
        static const T c3 = ( T ) 10./ ( T ) 3.;

        return
            ( val ( x+1,y ) +val ( x-1,y )
              + val ( x,y+1 ) +val ( x,y-1 ) ) *c1
            + ( val ( x+1,y+1 ) +val ( x-1,y-1 )
                +val ( x+1,y-1 ) +val ( x-1,y+1 ) ) *c2
            -val ( x,y ) *c3;
    }

private:
    ScalarField2D<T> const& field;
    plint x,y;
};


/* *************** Compute_C_processor ***************** */

template<typename T, template<typename U> class Descriptor >
Compute_C_processor2D <T,Descriptor>::Compute_C_processor2D ( T M_ )
    : M ( M_ )
{ }

template<typename T, template<typename U> class Descriptor >
void Compute_C_processor2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    BlockLattice2D<T,Descriptor>& f = *dynamic_cast<BlockLattice2D<T,Descriptor>*> ( blocks[0] );
    // Laplacian of chemical potential at previous time step t-1.
    ScalarField2D<T>& laplaceMu     = *dynamic_cast<ScalarField2D<T>*> ( blocks[1] );
    ScalarField2D<T>& C             = *dynamic_cast<ScalarField2D<T>*> ( blocks[2] );
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            // Use templates to compute order-0 moment of f.
            C.get ( iX,iY ) = momentTemplates<T,Descriptor>::get_rhoBar ( f.get ( iX,iY ) )
                              + 0.5 * M * laplaceMu.get ( iX,iY );
        }
    }
}


template<typename T, template<typename U> class Descriptor >
Compute_C_processor2D<T,Descriptor>*
Compute_C_processor2D<T,Descriptor>::clone() const
{
    return new Compute_C_processor2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor >
void Compute_C_processor2D<T,Descriptor>::getTypeOfModification (
    std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing; // f
    modified[1] = modif::nothing; // laplaceMu(t-1)
    modified[2] = modif::staticVariables;  // C
}


/* *************** Compute_gradC_rho_mu_processor ***************** */

template<typename T>
Compute_gradC_rho_mu_processor2D<T>::Compute_gradC_rho_mu_processor2D (
    T beta_, T kappa_, T rho_h_, T rho_l_ )
    : beta ( beta_ ),
      kappa ( kappa_ ),
      rho_h ( rho_h_ ),
      rho_l ( rho_l_ )
{ }

template<typename T>
void Compute_gradC_rho_mu_processor2D<T>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    // blocks[0] stands for the unused f.
    ScalarField2D<T>& C       = *dynamic_cast<ScalarField2D<T>*> ( blocks[1] );
    TensorField2D<T,2>& gradC = *dynamic_cast<TensorField2D<T,2>*> ( blocks[2] );
    ScalarField2D<T>& rho     = *dynamic_cast<ScalarField2D<T>*> ( blocks[3] );
    ScalarField2D<T>& mu      = *dynamic_cast<ScalarField2D<T>*> ( blocks[4] );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            FiniteDifference2D<T> ( C,iX,iY ).centralGradient ( gradC.get ( iX,iY ) );
            T C_ = C.get ( iX,iY );
            rho.get ( iX,iY ) = C_ * ( rho_h-rho_l ) + rho_l;
            mu.get ( iX,iY ) = 4.*beta*C_* ( C_-0.5 ) * ( C_-1. ) -
                               kappa*FiniteDifference2D<T> ( C,iX,iY ).laplacian();
        }
    }
}


template<typename T>
Compute_gradC_rho_mu_processor2D<T>*
Compute_gradC_rho_mu_processor2D<T>::clone() const
{
    return new Compute_gradC_rho_mu_processor2D<T> ( *this );
}

template<typename T>
void Compute_gradC_rho_mu_processor2D<T>::getTypeOfModification (
    std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;  // f
    modified[1] = modif::nothing;  // C
    modified[2] = modif::staticVariables;   // gradC
    modified[3] = modif::staticVariables;   // rho
    modified[4] = modif::staticVariables;   // mu
}


/* *************** Compute_gradMu_laplaceMu_processor ***************** */

template<typename T, template<typename U> class Descriptor >
Compute_gradMu_laplaceMu_processor2D <T,Descriptor>::Compute_gradMu_laplaceMu_processor2D()
{ }

template<typename T, template<typename U> class Descriptor >
void Compute_gradMu_laplaceMu_processor2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    // blocks[0] stands for the unused f.
    ScalarField2D<T>& mu            = *dynamic_cast<ScalarField2D<T>*> ( blocks[1] );
    TensorField2D<T,2>& gradMu      = *dynamic_cast<TensorField2D<T,2>*> ( blocks[2] );
    ScalarField2D<T>& laplaceMu     = *dynamic_cast<ScalarField2D<T>*> ( blocks[3] );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            FiniteDifference2D<T> ( mu,iX,iY ).centralGradient ( gradMu.get ( iX,iY ) );
            laplaceMu.get ( iX,iY ) = FiniteDifference2D<T> ( mu,iX,iY ).laplacian();
        }
    }
}


template<typename T, template<typename U> class Descriptor >
Compute_gradMu_laplaceMu_processor2D<T,Descriptor>*
Compute_gradMu_laplaceMu_processor2D<T,Descriptor>::clone() const
{
    return new Compute_gradMu_laplaceMu_processor2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor >
void Compute_gradMu_laplaceMu_processor2D<T,Descriptor>::getTypeOfModification (
    std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;  // f
    modified[1] = modif::nothing;  // mu
    modified[2] = modif::staticVariables;   // gradMu
    modified[3] = modif::staticVariables;   // laplaceMu
}


/* *************** Compute_gradMu_laplaceMu_u_p1_processor ***************** */

template<typename T, template<typename U> class Descriptor >
Compute_gradMu_laplaceMu_u_p1_processor2D <T,Descriptor>::Compute_gradMu_laplaceMu_u_p1_processor2D (
    T rho_h_, T rho_l_, T RT_ )
    : rho_h ( rho_h_ ),
      rho_l ( rho_l_ ),
      RT ( RT_ )
{ }

template<typename T, template<typename U> class Descriptor >
void Compute_gradMu_laplaceMu_u_p1_processor2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    // blocks[0] stands for the unused f.
    BlockLattice2D<T,Descriptor>& g = *dynamic_cast<BlockLattice2D<T,Descriptor>*> ( blocks[1] );
    ScalarField2D<T>& C             = *dynamic_cast<ScalarField2D<T>*> ( blocks[2] );
    ScalarField2D<T>& rho           = *dynamic_cast<ScalarField2D<T>*> ( blocks[3] );
    TensorField2D<T,2>& gradC       = *dynamic_cast<TensorField2D<T,2>*> ( blocks[4] );
    ScalarField2D<T>& mu            = *dynamic_cast<ScalarField2D<T>*> ( blocks[5] );
    TensorField2D<T,2>& gradMu      = *dynamic_cast<TensorField2D<T,2>*> ( blocks[6] );
    ScalarField2D<T>& laplaceMu     = *dynamic_cast<ScalarField2D<T>*> ( blocks[7] );
    TensorField2D<T,2>& u           = *dynamic_cast<TensorField2D<T,2>*> ( blocks[8] );
    ScalarField2D<T>& p1            = *dynamic_cast<ScalarField2D<T>*> ( blocks[9] );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Array<T,2>& u_      = u.get ( iX,iY );
            Array<T,2>& gradMu_ = gradMu.get ( iX,iY );
            Array<T,2>& gradC_  = gradC.get ( iX,iY );
            T invRho            = ( T ) 1 / rho.get ( iX,iY );

            FiniteDifference2D<T> ( mu,iX,iY ).centralGradient ( gradMu_ );
            laplaceMu.get ( iX,iY ) = FiniteDifference2D<T> ( mu,iX,iY ).laplacian();
            // Use templates to compute order-0 and order-1 moment of g.
            momentTemplates<T,Descriptor>::get_rhoBar_j (
                g.get ( iX,iY ), p1.get ( iX,iY ), u_ );
            u_ *= invRho/RT;
            u_ -= gradMu_ * ( 0.5*invRho*C.get ( iX,iY ) );
            p1.get ( iX,iY ) +=
                0.5*RT * ( rho_h-rho_l )
                * VectorTemplateImpl<T,2>::scalarProduct ( u_, gradC_ );
        }
    }
}


template<typename T, template<typename U> class Descriptor >
Compute_gradMu_laplaceMu_u_p1_processor2D<T,Descriptor>*
Compute_gradMu_laplaceMu_u_p1_processor2D<T,Descriptor>::clone() const
{
    return new Compute_gradMu_laplaceMu_u_p1_processor2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor >
void Compute_gradMu_laplaceMu_u_p1_processor2D<T,Descriptor>::getTypeOfModification (
    std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;  // f
    modified[1] = modif::nothing;  // g
    modified[2] = modif::nothing;  // C
    modified[3] = modif::nothing;  // rho
    modified[4] = modif::nothing;  // gradC
    modified[5] = modif::nothing;  // mu
    modified[6] = modif::staticVariables;   // gradMu
    modified[7] = modif::staticVariables;   // laplaceMu
    modified[8] = modif::staticVariables;   // u
    modified[9] = modif::staticVariables;   // p1
}


/* *************** HeLeeCollisionProcessor ***************** */

template<typename T, template<typename U> class Descriptor >
HeLeeCollisionProcessor2D <T,Descriptor>::HeLeeCollisionProcessor2D (
    T rho_h_, T rho_l_, T tau_h_, T tau_l_, T M_, T RT_,
    bool initialize_ )
    : rho_h ( rho_h_ ),
      rho_l ( rho_l_ ),
      tau_h ( tau_h_ ),
      tau_l ( tau_l_ ),
      M ( M_ ),
      RT ( RT_ ),
      initialize ( initialize_ )
{ }

template<typename T, template<typename U> class Descriptor >
void HeLeeCollisionProcessor2D<T,Descriptor>::computeAdvectionTerms (
    ScalarField2D<T> const& C, T& adv_gradC, T& bias_adv_gradC,
    plint iX, plint iY, plint iPop )
{
    typedef Descriptor<T> D;
    T diff_p2 = C.get ( iX+2*D::c[iPop][0], iY+2*D::c[iPop][1] ) -
                C.get ( iX+D::c[iPop][0], iY+D::c[iPop][1] );
    T diff_p1 = C.get ( iX+D::c[iPop][0], iY+D::c[iPop][1] ) -
                C.get ( iX, iY );
    T diff_p0 = C.get ( iX, iY ) -
                C.get ( iX-D::c[iPop][0], iY-D::c[iPop][1] );
    adv_gradC = diff_p1 - 0.5* ( diff_p1-diff_p0 );
    bias_adv_gradC = diff_p1 - 0.25* ( diff_p2-diff_p0 );
}

template<typename T, template<typename U> class Descriptor >
void HeLeeCollisionProcessor2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    typedef Descriptor<T> D;
    BlockLattice2D<T,Descriptor>& f = *dynamic_cast<BlockLattice2D<T,Descriptor>*> ( blocks[0] );
    BlockLattice2D<T,Descriptor>& g = *dynamic_cast<BlockLattice2D<T,Descriptor>*> ( blocks[1] );
    ScalarField2D<T>& C             = *dynamic_cast<ScalarField2D<T>*> ( blocks[2] );
    ScalarField2D<T>& rho           = *dynamic_cast<ScalarField2D<T>*> ( blocks[3] );
    TensorField2D<T,2>& gradC       = *dynamic_cast<TensorField2D<T,2>*> ( blocks[4] );
    ScalarField2D<T>& mu            = *dynamic_cast<ScalarField2D<T>*> ( blocks[5] );
    TensorField2D<T,2>& gradMu      = *dynamic_cast<TensorField2D<T,2>*> ( blocks[6] );
    ScalarField2D<T>& laplaceMu     = *dynamic_cast<ScalarField2D<T>*> ( blocks[7] );
    TensorField2D<T,2>& u           = *dynamic_cast<TensorField2D<T,2>*> ( blocks[8] );
    ScalarField2D<T>& p1            = *dynamic_cast<ScalarField2D<T>*> ( blocks[9] );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Array<T,2>& u_ = u.get ( iX,iY );
            T& C_ = C.get ( iX,iY );
            T& rho_ = rho.get ( iX,iY );
            T& laplaceMu_ = laplaceMu.get ( iX,iY );
            T& p1_ = p1.get ( iX,iY );
            Cell<T,Descriptor>& f_ = f.get ( iX,iY );
            Cell<T,Descriptor>& g_ = g.get ( iX,iY );

            T tau = C_* ( tau_h-tau_l ) + tau_l;

            Array<T,2> biasGradC, biasGradMu, gradP, biasGradP;
            FiniteDifference2D<T> ( C,iX,iY ).biasedGradient ( biasGradC );
            FiniteDifference2D<T> ( mu,iX,iY ).biasedGradient ( biasGradMu );
            FiniteDifference2D<T> ( p1,iX,iY ).centralGradient ( gradP );
            FiniteDifference2D<T> ( p1,iX,iY ).biasedGradient ( biasGradP );

            T uGradC =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, gradC.get ( iX,iY ) );
            T uBiasGradC =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, biasGradC );
            T uGradMu =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, gradMu.get ( iX,iY ) );
            T uBiasGradMu =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, biasGradMu );
            T uGradP =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, gradP );
            T uBiasGradP =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, biasGradP );
            T uSqr =
                VectorTemplateImpl<T,2>::normSqr ( u_ );
            for ( plint iPop=0; iPop<Descriptor<T>::q; ++iPop )
            {
                T ci_u = D::c[iPop][0]*u_[0]+D::c[iPop][1]*u_[1];
                T ti = Descriptor<T>::t[iPop];
                T gamma0 = ti;
                T gammaBar = 3.*ci_u + 4.5*ci_u*ci_u - 1.5*uSqr;
                T gamma = ti* ( 1.+gammaBar );

                T adv_gradC, bias_adv_gradC;
                computeAdvectionTerms ( C, adv_gradC, bias_adv_gradC, iX,iY,iPop );
                adv_gradC -= uGradC;
                bias_adv_gradC -= uBiasGradC;

                T adv_gradP, bias_adv_gradP;
                computeAdvectionTerms ( p1, adv_gradP, bias_adv_gradP, iX,iY,iPop );
                adv_gradP -= uGradP;
                bias_adv_gradP -= uBiasGradP;

                T adv_gradMu, bias_adv_gradMu;
                computeAdvectionTerms ( mu, adv_gradMu, bias_adv_gradMu, iX,iY,iPop );
                adv_gradMu -= uGradMu;
                adv_gradMu *= C_;
                bias_adv_gradMu -= uBiasGradMu;
                bias_adv_gradMu *= C_;

                T fieq = ti * C_* ( 1+gammaBar )
                         - 0.5*gamma* ( adv_gradC-1./RT* ( adv_gradP+adv_gradMu ) *C_/rho_ )
                         - 0.5*gamma*M*laplaceMu_;

                T gieq = ti* ( p1_+rho_*RT*gammaBar )
                         - 0.5*RT*adv_gradC* ( rho_h-rho_l ) * ( gamma-gamma0 )
                         + 0.5*gamma*adv_gradMu;
                if ( initialize )
                {
                    f_[iPop] = fieq;
                    g_[iPop] = gieq;
                }
                else
                {
                    f_[iPop] += - ( f_[iPop]-fieq ) * 1./ ( tau+0.5 )
                                +gamma* (
                                    bias_adv_gradC- ( bias_adv_gradP+bias_adv_gradMu ) *1./RT*C_/rho_ )
                                + M*laplaceMu_*gamma;
                    g_[iPop] += - ( g_[iPop]-gieq ) * 1./ ( tau+0.5 )
                                + bias_adv_gradC* ( rho_h-rho_l ) *RT* ( gamma-gamma0 )-bias_adv_gradMu*gamma;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor >
HeLeeCollisionProcessor2D<T,Descriptor>*
HeLeeCollisionProcessor2D<T,Descriptor>::clone() const
{
    return new HeLeeCollisionProcessor2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor >
void HeLeeCollisionProcessor2D<T,Descriptor>::getTypeOfModification (
    std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;   // f
    modified[1] = modif::staticVariables;   // g
    modified[2] = modif::nothing;  // C
    modified[3] = modif::nothing;  // rho
    modified[4] = modif::nothing;  // gradC
    modified[5] = modif::nothing;  // mu
    modified[6] = modif::nothing;  // gradMu
    modified[7] = modif::nothing;  // laplaceMu
    modified[8] = modif::nothing;  // u
    modified[9] = modif::nothing;  // p1
}

/* ****** HeLeeCollisionProcessor for ForcedD2Q9Descriptor ******** */

template<typename T >
HeLeeCollisionProcessor2D <T,ForcedD2Q9Descriptor>::HeLeeCollisionProcessor2D (
    T rho_h_, T rho_l_, T tau_h_, T tau_l_, T M_, T RT_,
    bool initialize_ )
    : rho_h ( rho_h_ ),
      rho_l ( rho_l_ ),
      tau_h ( tau_h_ ),
      tau_l ( tau_l_ ),
      M ( M_ ),
      RT ( RT_ ),
      initialize ( initialize_ )
{ }

template<typename T>
void HeLeeCollisionProcessor2D<T,ForcedD2Q9Descriptor>::computeAdvectionTerms (
    ScalarField2D<T> const& C, T& adv_gradC, T& bias_adv_gradC,
    plint iX, plint iY, plint iPop )
{
    typedef ForcedD2Q9Descriptor<T> D;
    T diff_p2 = C.get ( iX+2*D::c[iPop][0], iY+2*D::c[iPop][1] ) -
                C.get ( iX+D::c[iPop][0], iY+D::c[iPop][1] );
    T diff_p1 = C.get ( iX+D::c[iPop][0], iY+D::c[iPop][1] ) -
                C.get ( iX, iY );
    T diff_p0 = C.get ( iX, iY ) -
                C.get ( iX-D::c[iPop][0], iY-D::c[iPop][1] );
    adv_gradC = diff_p1 - 0.5* ( diff_p1-diff_p0 );
    bias_adv_gradC = diff_p1 - 0.25* ( diff_p2-diff_p0 );
}

template <typename T>
void HeLeeCollisionProcessor2D<T,ForcedD2Q9Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    typedef descriptors::ForcedD2Q9Descriptor<T> D;
    BlockLattice2D<T,ForcedD2Q9Descriptor>& f = *dynamic_cast<BlockLattice2D<T,ForcedD2Q9Descriptor>*> ( blocks[0] );
    BlockLattice2D<T,ForcedD2Q9Descriptor>& g = *dynamic_cast<BlockLattice2D<T,ForcedD2Q9Descriptor>*> ( blocks[1] );
    ScalarField2D<T>& C             = *dynamic_cast<ScalarField2D<T>*> ( blocks[2] );
    ScalarField2D<T>& rho           = *dynamic_cast<ScalarField2D<T>*> ( blocks[3] );
    TensorField2D<T,2>& gradC       = *dynamic_cast<TensorField2D<T,2>*> ( blocks[4] );
    ScalarField2D<T>& mu            = *dynamic_cast<ScalarField2D<T>*> ( blocks[5] );
    TensorField2D<T,2>& gradMu      = *dynamic_cast<TensorField2D<T,2>*> ( blocks[6] );
    ScalarField2D<T>& laplaceMu     = *dynamic_cast<ScalarField2D<T>*> ( blocks[7] );
    TensorField2D<T,2>& u           = *dynamic_cast<TensorField2D<T,2>*> ( blocks[8] );
    ScalarField2D<T>& p1            = *dynamic_cast<ScalarField2D<T>*> ( blocks[9] );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Array<T,2>& u_ = u.get ( iX,iY );
            T& C_ = C.get ( iX,iY );
            T& rho_ = rho.get ( iX,iY );
            T& laplaceMu_ = laplaceMu.get ( iX,iY );
            T& p1_ = p1.get ( iX,iY );
            Cell<T,ForcedD2Q9Descriptor>& f_ = f.get ( iX,iY );
            Cell<T,ForcedD2Q9Descriptor>& g_ = g.get ( iX,iY );

            T tau = C_* ( tau_h-tau_l ) + tau_l;

            Array<T,2> biasGradC, biasGradMu, gradP, biasGradP;
            FiniteDifference2D<T> ( C,iX,iY ).biasedGradient ( biasGradC );
            FiniteDifference2D<T> ( mu,iX,iY ).biasedGradient ( biasGradMu );
            FiniteDifference2D<T> ( p1,iX,iY ).centralGradient ( gradP );
            FiniteDifference2D<T> ( p1,iX,iY ).biasedGradient ( biasGradP );

            T uGradC =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, gradC.get ( iX,iY ) );
            T uBiasGradC =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, biasGradC );
            T uGradMu =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, gradMu.get ( iX,iY ) );
            T uBiasGradMu =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, biasGradMu );
            T uGradP =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, gradP );
            T uBiasGradP =
                VectorTemplateImpl<T,2>::scalarProduct ( u_, biasGradP );
            T uSqr =
                VectorTemplateImpl<T,2>::normSqr ( u_ );
            for ( plint iPop=0; iPop<ForcedD2Q9Descriptor<T>::q; ++iPop )
            {
                T ci_u = D::c[iPop][0]*u_[0]+D::c[iPop][1]*u_[1];
                T ti = ForcedD2Q9Descriptor<T>::t[iPop];
                T gamma0 = ti;
                T gammaBar = 3.*ci_u + 4.5*ci_u*ci_u - 1.5*uSqr;
                T gamma = ti* ( 1.+gammaBar );

                T adv_gradC, bias_adv_gradC;
                computeAdvectionTerms ( C, adv_gradC, bias_adv_gradC, iX,iY,iPop );
                adv_gradC -= uGradC;
                bias_adv_gradC -= uBiasGradC;

                T adv_gradP, bias_adv_gradP;
                computeAdvectionTerms ( p1, adv_gradP, bias_adv_gradP, iX,iY,iPop );
                adv_gradP -= uGradP;
                bias_adv_gradP -= uBiasGradP;

                T adv_gradMu, bias_adv_gradMu;
                computeAdvectionTerms ( mu, adv_gradMu, bias_adv_gradMu, iX,iY,iPop );
                adv_gradMu -= uGradMu;
                adv_gradMu *= C_;
                bias_adv_gradMu -= uBiasGradMu;
                bias_adv_gradMu *= C_;

                T fieq = ti * C_* ( 1+gammaBar )
                         - 0.5*gamma* ( adv_gradC-1./RT* ( adv_gradP+adv_gradMu ) *C_/rho_ )
                         - 0.5*gamma*M*laplaceMu_;

                T gieq = ti* ( p1_+rho_*RT*gammaBar )
                         - 0.5*RT*adv_gradC* ( rho_h-rho_l ) * ( gamma-gamma0 )
                         + 0.5*gamma*adv_gradMu;
                if ( initialize )
                {
                    f_[iPop] = fieq;
                    g_[iPop] = gieq;
                }
                else
                {
                    f_[iPop] += - ( f_[iPop]-fieq ) * 1./ ( tau+0.5 )
                                +gamma* (
                                    bias_adv_gradC- ( bias_adv_gradP+bias_adv_gradMu ) *1./RT*C_/rho_ )
                                + M*laplaceMu_*gamma;
                    g_[iPop] += - ( g_[iPop]-gieq ) * 1./ ( tau+0.5 )
                                + bias_adv_gradC* ( rho_h-rho_l ) *RT* ( gamma-gamma0 )-bias_adv_gradMu*gamma;
                }
            }
            /*
             * FIXME:不能使用Guo等的作用力处理模型，必须使用另外的处理方式
             */
            externalForceTemplates<T,ForcedD2Q9Descriptor>::addGuoForce ( f_, u_, f_.getDynamics().getOmega(), rho_ );
            externalForceTemplates<T,ForcedD2Q9Descriptor>::addGuoForce ( g_, u_, g_.getDynamics().getOmega(), rho_ );
        }
    }
}

template<typename T>
HeLeeCollisionProcessor2D<T,ForcedD2Q9Descriptor>*
HeLeeCollisionProcessor2D<T,ForcedD2Q9Descriptor>::clone() const
{
    return new HeLeeCollisionProcessor2D<T,ForcedD2Q9Descriptor> ( *this );
}

template<typename T>
void HeLeeCollisionProcessor2D<T,ForcedD2Q9Descriptor>::getTypeOfModification (
    std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;   // f
    modified[1] = modif::staticVariables;   // g
    modified[2] = modif::nothing;  // C
    modified[3] = modif::nothing;  // rho
    modified[4] = modif::nothing;  // gradC
    modified[5] = modif::nothing;  // mu
    modified[6] = modif::nothing;  // gradMu
    modified[7] = modif::nothing;  // laplaceMu
    modified[8] = modif::nothing;  // u
    modified[9] = modif::nothing;  // p1
}

}  // namespace plb

#endif  // HE_LEE_PROCESSOR_2D_HH
