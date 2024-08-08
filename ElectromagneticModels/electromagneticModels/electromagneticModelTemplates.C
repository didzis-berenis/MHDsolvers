/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

//#include "electromagneticModel.H"
//#include "wordIOList.H"
//#include "compileTemplate.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{

    //template<class electromagneticModel>
    Foam::autoPtr<electromagneticModel> electromagneticModel::New
    (
        const fvMesh& mesh,
        const word& phaseName
    )
    {
        const IOdictionary electromagneticDict
        (
            physicalProperties::findModelDict(mesh, phaseName)
        );
        const word modelType(electromagneticDict.lookup("electromagneticType"));

        Info<< "Selecting electromagnetics model " << modelType << endl;

        typename electromagneticModel::fvMeshConstructorTable::iterator
            cstrIter =
            electromagneticModel::fvMeshConstructorTablePtr_->find(modelType);

        if (cstrIter == electromagneticModel::fvMeshConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown " << electromagneticModel::typeName << " type "
                << modelType << nl << nl
                << "Valid " << electromagneticModel::typeName << " types are:" << nl
                << electromagneticModel::fvMeshConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return autoPtr<electromagneticModel>(cstrIter()(mesh, phaseName));
    }
}

// ************************************************************************* //
