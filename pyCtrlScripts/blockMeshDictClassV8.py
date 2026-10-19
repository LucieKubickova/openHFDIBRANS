#!/usr/bin/python

#FILE DESCRIPTION=======================================================

# Python class containing the functions to ease up generation of
# blockMeshDict files in OpenFOAM

#IMPORT BLOCK===========================================================
import math
import numpy as np

class mesh:
    """ python class with mesh data """

    def __init__(self, mScale = 1.0):
        self.mScale = mScale
        self.vertices = list()
        self.nPoints = 0
        self.blocks = list()
        self.patches = list()
        self.edges = list()

    def areSame(self, p1, p2):
        p1 = np.array(p1)
        p2 = np.array(p2)

        if np.linalg.norm(p1 - p2) < 1e-15:
            return True

        else:
            return False

    def addVertex(self, vertex, neighbours):
        duplicate = False

        if len(neighbours) > 0:
            for neighbour in neighbours:
                nIndices = neighbour.indices

                for nIndex in nIndices:
                    nVertex = self.vertices[nIndex]
                    duplicate = self.areSame(vertex, nVertex)

                    if duplicate:
                        break
                
                if duplicate:
                    break

        if not duplicate:
            for nIndex in range(self.nPoints):
                nVertex = self.vertices[nIndex]
                duplicate = self.areSame(vertex, nVertex)

                if duplicate:
                    break

        if duplicate:
            return nIndex

        else:
            self.vertices.append(vertex)
            self.nPoints += 1
            return self.nPoints-1

    def addBlock(self, vertices, neighbours, nCells, grading, grType = "simpleGrading", name = None):
        block = blockClass(neighbours, nCells, grading, grType, name)
        self.blocks.append(block)

        for vertex in vertices:
            index = self.addVertex(vertex, block.neighs)
            block.indices.append(index)

        return block

    def addPatch(self, patchName, patchType, faces, options = list()):
        patch = patchClass(patchName, patchType, faces, options = options)
        self.patches.append(patch)

    def addEdge(self, edgeType, indices, vertices):
        edge = edgeClass(edgeType, indices, vertices)
        self.edges.append(edge)

    def writeBMD(self, path):
        bMD = open(path + "/blockMeshDict",'w')                    #open file for writing
        
        #-------------------------------------------------------------------
        # write the headline
        bMD.write("/*--------------------------------*- C++ -*----------------------------------*\ \n")
        bMD.write("| =========                 |                                                 | \n")
        bMD.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n")
        bMD.write("|  \\    /   O peration     | Version:  4.1                                   | \n")
        bMD.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n")
        bMD.write("|    \\/     M anipulation  |                                                 | \n")
        bMD.write("\*---------------------------------------------------------------------------*/ \n")
        
        # write file description
        bMD.write("FoamFile \n")
        bMD.write("{ \n \t version \t 2.0; \n \t format \t ascii; \n")
        bMD.write(" \t class \t\t dictionary; \n \t object \t blockMeshDict; \n} \n")
        bMD.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n\n")
        
        #-------------------------------------------------------------------
        # convert to metres
        bMD.write("convertToMeters \t" + repr(self.mScale) + "; \n\n")
        
        #-------------------------------------------------------------------
        # write vertices
        bMD.write("vertices \n( \n")
        k = 0 
        for vertex in self.vertices:
            bMD.write("\t ( " + " ".join(str(e) for e in vertex) + " )\t//" + " %03d"%k + "\n")
            k = k+1 
        bMD.write("); \n\n")
        
        #-----------------------------------------------------------------------
        # write edges
        bMD.write("edges \n( \n")
        for edge in self.edges:
            for line in edge.retEdgeString():
                bMD.write(line)
        bMD.write("); \n\n")
        
        #-------------------------------------------------------------------
        # write blocks
        bMD.write("blocks \n( \n")
        for block in self.blocks:
            for line in block.retBlockString():
                bMD.write(line)
        bMD.write("); \n\n")
        
        #-------------------------------------------------------------------
        # write boundaries
        bMD.write("boundary \n( \n")
        for patch in self.patches:
            for line in patch.retBoundString():
                bMD.write(line)
        bMD.write("); \n\n")
        
        #-------------------------------------------------------------------
        # write empty mergePatchPairs
        bMD.write("mergePatchPairs \n( \n")
        bMD.write("); \n\n")
        
        #-------------------------------------------------------------------
        # close file
        bMD.close()

    def writeStitching(self, path, stitchPairs):
        """ function to automatically write script for mesh stitching """

        stitchSc = open(path + "/stitchMeshSc.sh",'w')                                       #open file for writing
    
        stitchSc.write('#!/bin/sh\n\n')
    
        k = 0
        for pair in stitchPairs:
            stitchSc.write("stitchMesh -perfect -overwrite " + pair[0] + " " + pair[1] + " >> log.stitchMesh_%d\n"%k)
            k +=1
    
        stitchSc.close()

class blockClass:
    """ Base python class containing functions to generate blockMeshDict
        blocks """
        
    def __init__(self, neighbours, nCells, grading, grType, name):
        self.neighs = neighbours
        self.nCells = nCells
        self.grading= grading
        self.grStr  = grType
        self.name   = name

        self.indices = list()
        
    # functions to return block faces
    def retFXY0(self):
        """ return vert. indices of face in XY plane with 
            lower Z coordinate"""
        ind = self.indices
        return [ind[0], ind[1], ind[2], ind[3]]
    
    def retFXYE(self):
        """ return vert. indices of face in XY plane with 
            higher Z coordinate"""
        ind = self.indices
        return [ind[4], ind[5], ind[6], ind[7]]
        
    def retFXZ0(self):
        """ return vert. indices of face in XZ plane with 
            lower Y coordinate"""
        ind = self.indices
        return [ind[0], ind[1], ind[5], ind[4]]
    
    def retFXZE(self):
        """ return vert. indices of face in XZ plane with 
            higher Y coordinate"""
        ind = self.indices
        return [ind[3], ind[2], ind[6], ind[7]]
        
    def retFYZ0(self):
        """ return vert. indices of face in YZ plane with 
            lower X coordinate"""
        ind = self.indices
        return [ind[0], ind[3], ind[7], ind[4]]
    
    def retFYZE(self):
        """ return vert. indices of face in YZ plane with 
            higher X coordinate"""
        ind = self.indices
        return [ind[1], ind[2], ind[6], ind[5]]

    # functions to return block edges
    def retEX0Y0(self):
        """ return vert. indices of an edge along Z axis
            with lower X coordinate and lower Y coordinate """
        ind = self.indices
        return [ind[0], ind[4]]
        
    def retEX0YE(self):
        """ return vert. indices of an edge along Z axis
            with lower X coordinate and higher Y coordinate """
        ind = self.indices
        return [ind[3], ind[7]]
        
    def retEXEYE(self):
        """ return vert. indices of an edge along Z axis
            with higher X coordinate and higher Y coordinate """
        ind = self.indices
        return [ind[2], ind[6]]
        
    def retEXEY0(self):
        """ return vert. indices of an edge along Z axis
            with higher X coordinate and lower Y coordinate """
        ind = self.indices
        return [ind[1], ind[5]]

    def retEX0Z0(self):
        """ return vert. indices of an edge along Y axis
            with lower X coordinate and lower Z coordinate """
        ind = self.indices
        return [ind[0], ind[3]]
        
    def retEX0ZE(self):
        """ return vert. indices of an edge along Y axis
            with lower X coordinate and higher Z coordinate """
        ind = self.indices
        return [ind[4], ind[7]]
        
    def retEXEZE(self):
        """ return vert. indices of an edge along Y axis
            with higher X coordinate and higher Z coordinate """
        ind = self.indices
        return [ind[5], ind[6]]
        
    def retEXEZ0(self):
        """ return vert. indices of an edge along Y axis
            with higher X coordinate and lower Z coordinate """
        ind = self.indices
        return [ind[1], ind[2]]

    def retEY0Z0(self):
        """ return vert. indices of an edge along X axis
            with lower Y coordinate and lower Z coordinate """
        ind = self.indices
        return [ind[0], ind[1]]
        
    def retEY0ZE(self):
        """ return vert. indices of an edge along X axis
            with lower Y coordinate and higher Z coordinate """
        ind = self.indices
        return [ind[4], ind[5]]
        
    def retEYEZE(self):
        """ return vert. indices of an edge along X axis
            with higher Y coordinate and higher Z coordinate """
        ind = self.indices
        return [ind[7], ind[6]]
        
    def retEYEZ0(self):
        """ return vert. indices of an edge along X axis
            with higher Y coordinate and lower Z coordinate """
        ind = self.indices
        return [ind[3], ind[2]]

    # functions to return the block for blockMeshDict
    def retBlockString(self):
        introStr    = "\thex\n"
        vertStr     = "\t\t(" + " ".join(" %d"%vert for vert in self.indices) + ")"
        if not self.name == None:
            nameStr     = "\t" + self.name
        discStr     = "\t(" + " ".join(" %d"%nC for nC in self.nCells) + ")"
        grTpStr     = "\t" + self.grStr + " "
        gradStr     = "\t(" + " ".join(" %s"%gr for gr in self.grading) + ")\n\n"

        if not self.name == None:
            return [introStr,vertStr,nameStr,discStr,grTpStr,gradStr]
        else:
            return [introStr,vertStr,discStr,grTpStr,gradStr]

class edgeClass:
    def __init__(self, edgeType, indices, vertices):
        self.edgeType = edgeType
        self.indices = indices
        self.vertices = vertices

    def retEdgeString(self):
        typeStr     = "\t"+ self.edgeType + " "
        indStr      = str(self.indices[0]) + " " + str(self.indices[1]) + "\n\t"

        if self.edgeType == "polyLine":
            indStr += "\t(\n"

        vertStr = ""
        for vertex in self.vertices:
            vertStr += "\t( " + " ".join(str(e) for e in vertex) + " )\n"

        if self.edgeType == "polyLine":
            vertStr += "\t)\n"

        outStr = [typeStr, indStr, vertStr]

        return outStr

class patchClass:
    def __init__(self, patchName, patchType, faces, options = list()):
        self.patchName = patchName
        self.patchType = patchType
        self.faces = faces
        self.options = options

    def retBoundString(self):
        nameStr     = "\t" + self.patchName + "\n\t{\n"
        typeStr     = "\t\ttype " + self.patchType + ";\n"
        optStr      = ""
        for option in self.options:
            optStr += "\t\t" + option[0] + " " + option[1] + ";\n"
        faceStr0    = "\t\tfaces\n\t\t(\n"
        faceStrE    = "\t\t);\n\t}\n\n"
        outStr      = [nameStr,typeStr,optStr,faceStr0]

        for face in self.faces:
            outStr.append("\t\t\t(" + " ".join(" %d"%fc for fc in face) + ")\n")

        outStr.append(faceStrE)

        return outStr
