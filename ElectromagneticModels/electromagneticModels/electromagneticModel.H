/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::electromagneticModel

Description
    Base-class for fluid and solid thermodynamic properties.

    The electromagneticModel class is the pure virtual interface, plus static data and
    functions, including the selection table. It is default-constructable and
    holds no non-static data.

    The electromagneticModel::implementation class contains the class data, non-default
    construction and other implementation details.

    This structure allows for arbitrary mixing of the interfaces. Derived
    levels should implement a similar separation of interface and
    implementation. All interface classes should be default-constructable, and
    should be inherited virtually into whatever interfaces are required.
    Implementation classes should virtually inherit from their corresponding
    interface class and *not* a lower level implementation class.

    In derived levels that are complete enough to act as the base of a
    thermodynamic instantiation, an additional *::composite sub class should be
    defined, which (non-virtually) inherits from all the *::implementation
    classes that it needs. This way, the composite class composes an interface
    for which diamond patterns are resolved by virtual inheritance, whilst also
    obtaining a single copy of all the implementation classes it needs to form
    a complete implementation. The use of virtual inheritance does not result
    in additional constructor calls propagating further down the hierarchy
    (into heThermo and similar) because all virtually inherited interface
    classes are default constructable.

SourceFiles
    electromagneticModel.C

\*---------------------------------------------------------------------------*/

#ifndef electromagneticModel_H
#define electromagneticModel_H

#include "volFields.H"
#include "typeInfo.H"
#include "physicalProperties.H"
#include "autoPtr.H"
#include "wordIOList.H"
//#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class electromagneticModel Declaration
\*---------------------------------------------------------------------------*/

class electromagneticModel
:
    public physicalProperties
{
    // Private Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        // Switch for updating electromagnetics

        bool correctElectromagnetics_ = false;

protected:

    // Protected data

        //- Phase-name
        const word& phaseName_;


        // Fields

        //- Electrical conductivity [S/m]
        volScalarField sigma_;
        dimensionedScalar sigmaConst_;

        // Lorentz force

        volVectorField JxB_;

        // Joule heating

        volScalarField JJsigma_;

        //- Inverse Electric conductivity [m/S]
        volScalarField sigmaInv_;

    // Protected Static Member Functions

        //- Lookup and the named field, or construct it as MUST-READ if it is
        //  not found, and return a reference
        static volScalarField& lookupOrConstructScalar
        (
            const fvMesh& mesh,
            const char* name
        );

        //- Lookup and the named field, or construct it as MUST-READ if it is
        //  not found, and return a reference
        static volVectorField& lookupOrConstructVector
        (
            const fvMesh& mesh,
            const char* name
        );

        //- Construct named field as MUST-READ if it is not found
        static void constructScalar
        (
            const fvMesh& mesh,
            const char* name
        );

        //- Construct named field as MUST-READ if it is not found
        static void constructVector
        (
            const fvMesh& mesh,
            const char* name
        );

        //- Get the constructor iterator for the given thermo dictionary and
        //  constructor iterator table
        /*
        template<class ElectromagneticModel, class Table>
        static typename Table::iterator lookupCstrIter
        (
            const dictionary& thermoDict,
            Table* tablePtr
        );
        */
        //- Specialisation of the Generic New for electromagneticModel
        /*static autoPtr<electromagneticModel> New
        (
            const fvMesh& mesh,
            const word& phaseName
        );*/


public:

    //- Runtime type information
    TypeName("electromagneticModel");


    //- Declare run-time constructor selection table
    declareRunTimeNewSelectionTable
    (
        autoPtr,
        electromagneticModel,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    //- Generic New for each of the related electromagnetics model
    //template<class ElectromagneticModel>
    static autoPtr<electromagneticModel> New
    (
        const fvMesh& mesh,
        const word& phaseName=word::null
    );

    // Static Member data
        //- Name of the thermophysical properties dictionary
        //using physicalProperties::typeName;

    // Static Member Functions

        //- Lookup the thermo associated with the given patch field
        static const electromagneticModel& lookupElectromagnetic(const fvPatchScalarField& pf);


    // Constructors

        //- Construct from components
        electromagneticModel
        (
            const fvMesh& mesh,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        electromagneticModel(const electromagneticModel&) = delete;


    //- Destructor
    virtual ~electromagneticModel()
    {}


    // Member Functions

        // Lorentz force
        volVectorField& JxB;

        // Joule heating
        volScalarField& JJsigma;

        //- Properties dictionary
        virtual const IOdictionary& properties() const = 0;

        //- Read model coefficients if they have changed
        virtual bool read() = 0;

        //- Non-const access the properties dictionary
        virtual IOdictionary& properties() = 0;

        //- Return const access to the mesh
        virtual const fvMesh& mesh() const = 0;

        //- Phase name
        virtual const word& phaseName() const = 0;

        //- Inverse Electric conductivity [m/S]
        virtual const volScalarField sigmaInv() const = 0;

        //- Electric conductivity [S/m]
        virtual volScalarField sigma() = 0;

        //- Electric conductivity [S/m]
        virtual const volScalarField sigma() const = 0;

        //- Electric conductivity for patch [S/m]
        virtual scalarField sigma(const label patchi) const;

        //- Electric conductivity [S/m]
        virtual const dimensionedScalar sigmaConst() const = 0;

        //- Electric potential [V]
        virtual volScalarField& PotE(bool imaginary=false) = 0;

        //- Electric potential [V]
        virtual const volScalarField& PotE(bool imaginary=false) const = 0;

        //- Current density [A/m^2]
        virtual volVectorField& J(bool imaginary=false) = 0;

        //- Current density [A/m^2]
        virtual const volVectorField& J(bool imaginary=false) const = 0;

        //- Magnetic field [T]
        virtual volVectorField& B(bool imaginary=false) = 0;

        //- Magnetic field [T]
        virtual const volVectorField& B(bool imaginary=false) const = 0;

        //- Set correctElectromagnetics switch to true
        virtual void setCorrectElectromagnetics();

        //- Set correctElectromagnetics switch to false
        virtual void setCorrected();

        //- Get correctElectromagnetics switch value
        virtual bool correctElectromagnetics() const;

        //- return field from object registry
        virtual volVectorField& getVectorFromRegistry(const char* name);

        //- Name of the thermo physics
        virtual word electromagneticName() const = 0;

        //- Predict the electromagnetic model coefficients if possible
        //  without solving electromagnetic model equations
        virtual void predict() = 0;

        //- Solve the electromagnetic model equations
        //- Update properties
        virtual void correct() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "electromagneticModelTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
