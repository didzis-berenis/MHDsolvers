/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::solvers::incompressibleConductingVoFEtransfer

Description
    Solver module incompressibleConductingVoF with an optional
    electric-field transfer from the harmonic field solver. Select it in
    system/controlDict with
    \verbatim
        solver incompressibleConductingVoFEtransfer;
    \endverbatim

    Received field (constant/electromagneticModel):
    \verbatim
        receivedField       J;     // J (default) | E
        eTransferProjection none;  // none (default) | full | incremental
    \endverbatim

    receivedField J, or the key absent: the fields received from the
    harmonic field solver at a coupling exchange are the current density,
    rescaled between exchanges by the alpha1 ratio. This module then runs the
    incompressibleConductingVoF code unchanged.

    receivedField E (harmonic model only): the received fields are the
    electric field E (real and imaginary part). The module forms the current
    from the CURRENT alpha1 and U at every time step, with E and B held
    between exchanges:
    \verbatim
        W     = E + (U - U_exchange) x B
        J     = sigma(alpha1) W                        (projection none)
        JxB   = 0.5 (J_re x B_re + J_im x B_im)
        Joule = 0.5 sigma(alpha1) (|W_re - grad phi_re|^2
                                 + |W_im - grad phi_im|^2)
    \endverbatim
    with sigma(alpha1) = clip(alpha1, 0, 1) sigma. The Joule heat is formed
    from the field, so it does not divide a received current by a small
    mixture conductivity near the interface; this is the reason to use the E
    option when heating matters.

    Motional-term convention: Elmer's MagnetoDynamicsCalcFields adds
    v x curl(A) to the exported electric field when a Lorentz velocity is
    given, so the received E already contains U_exchange x B. Only the
    velocity change since the exchange is added (U_exchange = the velocity
    stored at the exchange, the one sent to the field solver).

    eTransferProjection:
      - none (default): J = sigma(alpha1) W, phi = 0. This is the validated
        mode: it brought the module Joule heat to within a few percent of the
        exact value on an AC sphere benchmark (J transfer: +9 to +21 %) and
        to within 3 % of the field solver's own Joule power in dynamic
        levitated-drop runs, with the force balance unchanged.
      - full: divergence-free projection. A potential phi solves
        laplacian(sigmaS, phi) = div((sigmaS W)_f . S_f), sigmaS =
        (floor + (1 - floor) a) sigma, zero normal flux on non-coupled
        walls; J is reconstructed conservatively from the projected face
        flux. On the sphere benchmark it removed the field solver's
        divergence error in the conductor core but put radial current into
        the diffuse interface band; not recommended.
      - incremental: only the divergence created since the exchange is
        projected out (sigma E at the exchange). Exploratory.
    The full and incremental modes solve for eTransferPhiRe/eTransferPhiIm,
    which need an fvSolution entry, e.g.
    \verbatim
        "eTransferPhi.*" { solver PCG; preconditioner DIC; tolerance 1e-10;
                           relTol 0; maxIter 20000; }
    \endverbatim
    Every projection logs an "eTransfer:" line per part; phi is copied to
    PotEre/PotEim and the received field is written as eTransferERe/Im.

    Elmer side for receivedField E (case.sif): in the
    MagnetoDynamicsCalcFields solver set
    \verbatim
        Calculate Electric Field = Logical True
    \endverbatim
    (keep Calculate Current Density and Calculate Joule Heating on if the
    Elmer current or Joule output is used), and in the Elmer2OpenFOAM solver
    replace the six current-density targets by the electric field, in the
    same order:
    \verbatim
        Target Variable 1 = Electric Field re 1
        Target Variable 2 = Electric Field re 2
        Target Variable 3 = Electric Field re 3
        Target Variable 4 = Electric Field im 1
        Target Variable 5 = Electric Field im 2
        Target Variable 6 = Electric Field im 3
    \endverbatim
    (Magnetic Flux Density re/im 1-3 unchanged). foamRunEpot receives the
    vectors in a fixed order into Jre, Jim, Bre, Bim, so E arrives in the
    containers that hold J with receivedField J; the solver binary is the
    same for both options. With receivedField E the static (DC) field files
    are not handled, and the per-step potential (PotE) solve is skipped.

    The module derives from incompressibleConductingVoF and overrides only
    solveElectromagnetics, fixCurrentsWithMixtureConductivity and
    storeReceivedCurrents; with receivedField J each override calls the
    incompressibleConductingVoF implementation and nothing else.

SourceFiles
    incompressibleConductingVoFEtransfer.C

See also
    Foam::solvers::incompressibleConductingVoF

\*---------------------------------------------------------------------------*/

#ifndef incompressibleConductingVoFEtransfer_H
#define incompressibleConductingVoFEtransfer_H

#include "incompressibleConductingVoF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
              Class incompressibleConductingVoFEtransfer Declaration
\*---------------------------------------------------------------------------*/

class incompressibleConductingVoFEtransfer
:
    public incompressibleConductingVoF
{
    // Private Member Functions

        //- The receivedField key of the electromagnetic model (default J)
        word receivedField() const;


public:

    //- Runtime type information
    TypeName("incompressibleConductingVoFEtransfer");


    // Constructors

        //- Construct from region mesh
        incompressibleConductingVoFEtransfer(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        incompressibleConductingVoFEtransfer
        (
            const incompressibleConductingVoFEtransfer&
        ) = delete;


    //- Destructor
    virtual ~incompressibleConductingVoFEtransfer();


    // Member Functions

        //- Solve electromagnetic modelling
        //  (receivedField E: form the current from the received E)
        virtual void solveElectromagnetics();

        //- Fix currents based on mixture conductivity
        //  (receivedField E: form the current from the received E, alpha1
        //  and U, once per time step and received field)
        virtual void fixCurrentsWithMixtureConductivity();

        //- Store the just-received current density (and alpha1) as the
        //  reference for fixCurrentsWithMixtureConductivity
        //  (receivedField E: also invalidates the cached current)
        virtual void storeReceivedCurrents(bool imaginary=false);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const incompressibleConductingVoFEtransfer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
