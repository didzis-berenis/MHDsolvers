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

\*---------------------------------------------------------------------------*/

#include "conductingRegionSolvers.H"
#include <set>
#include <vector>
#include <map>

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

// Construct new global mesh from region meshes
Foam::autoPtr<Foam::fvMesh> Foam::conductingRegionSolvers::globalMeshNew_()
{
    // Prepare global mesh
    // Get points from each region
	List<vector> list_p;// Don't know a priori how many unique points there are
    int allPointI = 0;
    //First label comes from iterating through names_/faces/thisFacePoints
    std::map<label,label> allToUniquePointId;
    forAll(names_, i)
    {
        fvMesh& regionMesh = regions_[i];
        const faceList faces = regionMesh.faces();
        const pointField points = regionMesh.points();
        forAll(faces, faceI)
        {
            pointField thisFacePoints = faces[faceI].points(points);
            forAll (thisFacePoints,pointI)
            {
                point thisPoint = thisFacePoints[pointI];
                bool found = false;
                forAll(list_p, p)
                {
                    if (list_p[p] == thisPoint)
                    {
                        allToUniquePointId[allPointI] = p;
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    allToUniquePointId[allPointI] = list_p.size();
                    list_p.append(thisPoint);
                }
                allPointI++;
            }
        }
    }
    // Get faces from each region
    int face_size = 0;
    forAll(names_, i)
    {
        fvMesh& regionMesh = regions_[i];
        const faceList faces = regionMesh.faces();
        face_size += regionMesh.faces().size();
    }
    std::vector<face> reordered_faces;
    reordered_faces.reserve(face_size);//size known
    std::map<face,label> faceOwners;
    std::map<std::pair<word, label>, face> regionIdToFace;
    allPointI = 0;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const faceList faces = regionMesh.faces();
        const labelList owner = regionMesh.faceOwner();//including boundary faces
        const pointField points = regionMesh.points();
        forAll(faces, faceI)
        {
            pointField thisFacePoints = faces[faceI].points(points);
            labelList thisFaceLabels;
            thisFaceLabels.setSize(thisFacePoints.size());
            forAll (thisFacePoints,pointI)
            {
                thisFaceLabels[pointI] = allToUniquePointId[allPointI];//Get global point Id
                allPointI++;
            }
            face newFace(thisFaceLabels);
            std::pair<word, label> faceMarker = std::make_pair(regionName,faceI);//Marker for mapping to Face
            regionIdToFace[faceMarker] = newFace;
            reordered_faces.push_back(newFace);
            label thisOwnerCellID = owner[faceI];//Owner id in region mesh
            std::pair<word, label> cellMarker = std::make_pair(regionName,thisOwnerCellID);//Marker for mapping to Cell
            faceOwners[newFace] = regionToGlobalCellId[cellMarker];
            
        }
    }
    // Check if all faces have owners
    if (faceOwners.size() != reordered_faces.size())
    {
        FatalIOError
        << "Incorrect face owner size!\n"
        << exit(FatalIOError);
    }
    // Assign owners to boundary faces and neighbours to internal faces
    std::map<face,label> faceNeighbours;
    forAll(names_, i)
    {
        const word& regionName = names_[i].first();
        fvMesh& regionMesh = regions_[i];
        const cellList cells = regionMesh.cells();
        forAll(cells, cellI)
        {
            std::pair<word, label> cellMarker = std::make_pair(regionName,cellI);//Marker for mapping to Cell
            label globalCellID = regionToGlobalCellId[cellMarker];//Cell id in global mesh
            cell thisCell = cells[cellI];
            forAll (thisCell,faceI)
            {
                label regionFaceID = thisCell[faceI];//Face id in region mesh
                std::pair<word, label> faceMarker = std::make_pair(regionName,regionFaceID);//Marker for mapping to Face
                face thisFace = regionIdToFace[faceMarker];//Mapping to Face
                if (faceOwners[thisFace] != globalCellID)//This cell is not an owner of this face
                {
                    faceNeighbours[thisFace] = globalCellID;//Asign neighbour to internal face
                }
            }
        }
    }
    // Sort faces
    // First sort by owner, then by neighbour
    std::sort(reordered_faces.begin(), reordered_faces.end(),
    [&faceOwners](face const& f1, face const& f2)
        {
            return faceOwners[f1] < faceOwners[f2];
        }
    );
    std::sort(reordered_faces.begin(), reordered_faces.end(),
    [&faceNeighbours](face const& f1, face const& f2)
        {
            if (faceNeighbours.count(f1) && faceNeighbours.count(f2))//both faces have neighbours assigned
            {
                return faceNeighbours[f1] < faceNeighbours[f2];
            }
            else if (faceNeighbours.count(f1) && !faceNeighbours.count(f2))
            {
                //Prioritize existing neighbours over non-existing neighbours
                // (i.e. internal faces over boundary faces)
                return true;
            }
            return false;
        }
    );
    // Assign global face, owner, neighbour to lists.
    // Owner size is equal to face size,
    // because each face has an owner.
    // Neighbour size is less than face size,
    // because boundary faces have no neighbours.
	faceList list_f;
    list_f.setSize(reordered_faces.size());
    labelList list_o;
    list_o.setSize(reordered_faces.size());
    labelList list_n;
    int globalFaceI = 0;
    for (const auto& thisFace : reordered_faces)
    {
        list_f[globalFaceI] = thisFace;
        list_o[globalFaceI++] = faceOwners[thisFace];
        //auto it = faceNeighbours.find(thisFace);
        //if (it != faceNeighbours.end())
        if (faceNeighbours.count(thisFace))//face has neighbour
        {
            list_n.append(faceNeighbours[thisFace]);
        }
    }
    // Create global mesh pointer
    autoPtr<fvMesh> globalMeshPtr(
        new fvMesh
        (
            IOobject
            (
                // Using default region name ("region0").
                // Use different name if conflicts arise.
                fvMesh::defaultRegion,//"dummyRegion",//
                runTime_.timeName(),
                runTime_,
                IOobject::NO_READ
            ),
            pointField(list_p),
            faceList(list_f),
            labelList(list_o),
            labelList(list_n),
            false// No parallel sync
        )
        );
    fvMesh& globalMesh = globalMeshPtr();
    // Add boundary faces.
    // Create one boundary patch with all faces.
    // Size of the patch is the difference.
    // between the size of the owner list and the size of the neighbour list.
    // This is because boundary faces have no neighbours.
    // Faces were sorted, so that the boundary faces ar at the end of the list
    // (i.e. start index is the size of the neighbour list).
    List<polyPatch*> patches(1);
    patches[0] = polyPatch::New
            (
                "patch",//const word& patchType,
                "boundaryPatch",//const word& name,
                list_o.size()-list_n.size(),//const label size,
                list_n.size(),//const label start,
                0,//const label index,
                globalMesh.boundaryMesh()//const polyBoundaryMesh& bm
            ).ptr();
    // Add patches; no parallel comms
    globalMesh.addFvPatches(patches, false);
    // Get fvSchemes dictionary and transfer ownership to objectRegistry
    if (!globalMesh.objectRegistry::foundObject<IOdictionary>("fvSchemes"))
    {
        IOdictionary* fPtr
        (
            new IOdictionary
            (
                
            IOobject
            (
                "fvSchemes",
                globalMesh.time().system(),
                globalMesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
            )
        );
        fPtr->store(fPtr);
    }
    // Write global mesh for debugging purposes
    //globalMesh.write();
    return globalMeshPtr;
}

// ************************************************************************* //
