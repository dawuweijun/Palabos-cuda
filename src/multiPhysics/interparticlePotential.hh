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

#ifndef INTERPARTICLE_POTENTIAL_HH
#define INTERPARTICLE_POTENTIAL_HH

#include "multiPhysics/interparticlePotential.h"

namespace plb {

namespace interparticlePotential {

template<typename T>
PsiFunction<T>::~PsiFunction()
{ }

template<typename T>
T PsiIsRho<T>::compute(T rho) const {
    return rho;
}
template<typename T>
PsiIsRho<T>* PsiIsRho<T>::clone() const {
    return new  PsiIsRho<T>(*this);
}

template<typename T>
PsiShanChen93<T>::PsiShanChen93(T rho_0_)
    : rho_0(rho_0_)
{ }

template<typename T>
T PsiShanChen93<T>::compute(T rho) const {
    return rho_0 * ((T)1 - exp(-rho/rho_0));
}

template<typename T>
PsiShanChen93<T>* PsiShanChen93<T>::clone() const {
    return new PsiShanChen93<T>(*this);
}

template<typename T>
PsiShanChen94<T>::PsiShanChen94(T psi_0_, T rho_0_)
    : psi_0(psi_0_),
      rho_0(rho_0_)
{ }

template<typename T>
T PsiShanChen94<T>::compute(T rho) const {
    return psi_0 * exp(-rho_0/rho);
}

template<typename T>
PsiShanChen94<T>* PsiShanChen94<T>::clone() const {
    return new PsiShanChen94<T>(*this);
}

template<typename T>
PsiQian95<T>::PsiQian95(T rho_0_, T g_)
    : rho_0(rho_0_),
      rho_0_sqr(rho_0*rho_0),
      g(g_)
{ }

template<typename T>
T PsiQian95<T>::compute(T rho) const {
    T rho_sqr = rho*rho;
    return g * rho_0_sqr * rho_sqr /
             ( (T)2*(rho_0_sqr+rho_sqr+(T)2*rho_0*rho) );
}
template<typename T>
PsiQian95<T>* PsiQian95<T>::clone() const {
    return new PsiQian95<T>(*this);
}

}  // namespace interparticlePotential

}  // namespace plb

#endif  // INTERPARTICLE_POTENTIAL_HH
