# -- Python class to store, change and write OpenFOAM case data

# --    usage:
# --        1) load OpenFOAMCase from baseCase or from the different OpenFOAMCase
# --            -- loadOFCaseFromBaseCase(baseDir)
# --                    baseDir -- base case directory
# --        2) (optional) change OpenFOAMCase.dir (where the changes will be made), default -- baseDir
# --            -- changeOFCaseDir(dir)
# --        3) (optional) copy OpenFOAMCase.basedir to OpenFOAMCase.dir
# --            -- copyBaseCase()
# --        4) specify what you want to modify in openFoam source files:
# --            options:            
# --                a) replace([[in inFl (file), [whats (string)], [bys (string)]],...])
# --                        inFl -- openFoam source file you want to change
# --                        whats -- list of what you want to replace (string)
# --                        bys -- list of new values (string)
# --                b) setParameters([[in inFl (file), par (string), val (string), inSubDict (string)],...])
# --                        inFl -- openFoam source file you want to change
# --                        par -- parameter to change (string)
# --                        val -- new value of the parameter (string)
# --                        inSubDict -- subdictionary of the dictionary (string), if not used give ''
# --                c) addToDictionary([[in inFl (file), strToAdd (string), inSubDict (string)],...])
# --                        inFl -- openFoam source file you want to change
# --                        strToAdd -- string to add to dictionary
# --                        inSubDict -- subdictionary of the dictionary (string), if not used give ''
# --        6) run commands in OpenFOAMCase.dir
# --            -- runCommands(commands)

# --    additional functions:
# --        -- updateTimes() 
# --                update OpenFOAMCase.times variable (includes float folders in OpenFOAMCase.dir)
# --                and OpenFOAMCase.latestTime variable (max of OpenFOAMCase.times)

# -- standard imports 
import numpy as np
import os
import shutil as sh
import re

# -- custom function
from myAddFcs import *

class OpenFOAMCase:
    def __init__(self):
        """initialization of the new OpenFOAMCase"""
        # -- save starting directory
        self.whereIStart = os.getcwd()
        print("New OpenFOAMCase has been initialized, starting directory is:\n\t%s"%self.whereIStart)
    
    def loadOFCaseFromBaseCase(self, baseDir):
        """load OpenFOAMCase from base case directory"""
        print("OpenFOAMCase.baseCase has been set to :\n\t%s" % baseDir)
        self.baseDir = baseDir
        # -- OpenFOAMCase directory is initialized as its baseDir, can be changed by changeOFCaseDir(dir)
        self.dir = baseDir
    
    def changeOFCaseDir(self,dir):
        """change OpenFOAMCase directory"""
        self.dir = dir
        print("OpenFOAMCase directory has been changed from:\n\t%s to %s"%(self.baseDir,dir))

    def copyBaseCase(self):
        """copy OpenFOAMCase.baseDir to OpenFOAMCase.dir"""
        # -- if OpenFOAMCase.dir exists, remove it
        if os.path.isdir(self.dir): 
            sh.rmtree(self.dir)

        # -- copy baseDir to dir
        sh.copytree(self.baseDir,self.dir)
        print("OpenFOAMCase.baseDir has been copied to OpenFOAMCase.dir:\n\t %s --> %s"%(self.baseDir,self.dir))

    def replace(self, replaces):
        """replace option -- replaces = [[in inFl (file), [whats (string)], [bys (string)]]]"""
        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        # -- make the replaces
        for replace in replaces:
            inFl, what, by = replace
            with open(inFl, 'r') as fl:
                linesInFl = fl.readlines()
            for lnI in range(len(linesInFl)):
                for whatI in range(len(what)):
                    whatTu = what[whatI]
                    byTu = by[whatI]
                    if whatTu in linesInFl[lnI]:
                        linesInFl[lnI] = linesInFl[lnI].replace(whatTu,byTu)
                        print("In %s, I have replaced %s by %s on line %d." % (inFl, whatTu, byTu, lnI))
            with open(inFl, 'w') as fl:
                for lnI in range(len(linesInFl)):
                    fl.writelines(linesInFl[lnI])
        
        # -- move back where I start
        os.chdir(self.whereIStart)
    
    def setParameters(self, inParVals):
        """setParameter option -- inParVals = [[in inFl (file), par (string), val (string), inSubDict (string)]]"""
        """inSubDict option enables to find parameter in subDictionary, if "" then does nothing """
        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        # -- make the replaces
        for inParVal in inParVals:
            inFl, par, val, inSubDict = inParVal
            print('In %s, trying to change parameter %s in subDict %s:' % (inFl, par, inSubDict))
            with open(inFl, 'r') as fl:
                content = fl.read()
                contentN = re.sub(r'(%s)(\s+%s\s+[^\]]*]\s+|\s+)([\S ]+;)' % (par,par), r'\1 \2 %s;' % (val), content)
                matches = re.finditer(r'(%s)(\s+%s\s+[^\]]*]\s+|\s+)([\S ]+;)' % (par,par), content)
            if not inSubDict == "":
                match1 = re.search(r'%s\s*{([^}^{]*)[}{]' % inSubDict, content)
                contentMatchN = re.sub(r'(%s)(\s+%s\s+[^\]]*]\s+|\s+)([\S ]+;)' % (par,par), r'\1 \2 %s;' % (val), match1.group(0))
                matches = re.finditer(r'(%s)(\s+%s\s+[^\]]*]\s+|\s+)([\S ]+;)' % (par,par), match1.group(0))
                contentN = re.sub(r'%s\s*{([^}^{]*)[}{]' % inSubDict, contentMatchN, content)
            for match in (matches):
                print('\tIn %s, I have change val of %s parameter from %s to %s in subDict %s.' % (inFl, par, match.group(3), val, inSubDict))
            with open(inFl, 'w') as fl:
                fl.write(contentN)

        # -- move back where I start
        os.chdir(self.whereIStart)
    
    def addToDictionary(self, inParVals, brack = '{}'):
        """addToDictionary option -- inParVals = [[in inFl (file), strToAdd (string), inSubDict (string)]]"""
        """inSubDict option enables to find parameter in subDictionary, if "" then does nothing """
        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)
        
        # -- make the replaces
        for inParVal in inParVals:
            inFl, strToAdd, inSubDict = inParVal
            with open(inFl, 'r') as fl:
                linesInFl = fl.readlines()
            
            newLinesInFl = []
            # -- if inSubDict speciefied find where it is
            if not inSubDict == "":
                subDictSt, subDictEnd = -1, -1
                for lnI in range(len(linesInFl)):
                    if subDictSt != -1 and brack[0] in linesInFl[lnI]:
                        subDictEnd += 1
                    if subDictSt != -1 and brack[1] in linesInFl[lnI]:
                        subDictEnd -= 1
                    if inSubDict in linesInFl[lnI] and not '(' in linesInFl[lnI]:
                        subDictSt = lnI
                    elif (subDictSt != -1 and subDictEnd == -1):
                        subDictEnd = lnI
                if subDictSt == -1:
                    print("I could not find subDict %s in file %s."%(inSubDict,inFl))
                
                for lnI in range(len(linesInFl)):
                    if lnI < subDictEnd:
                        newLinesInFl.append(linesInFl[lnI])
                    elif lnI == (subDictEnd):
                        newLinesInFl.append('\t%s' % strToAdd)
                        print('In %s, I have added line %s to subdictionary %s on line %d' % (inFl, strToAdd, inSubDict, lnI))
                        newLinesInFl.append(linesInFl[lnI])
                    else:
                        newLinesInFl.append(linesInFl[lnI])
            else:
                newLinesInFl = linesInFl
                newLinesInFl.append('\n%s' % strToAdd)
                print('In %s, I have added line %s at the end.' % (inFl, strToAdd))
            with open(inFl, 'w') as fl:
                for lnI in range(len(newLinesInFl)):
                    fl.writelines(newLinesInFl[lnI])
        
        # -- move back where I start
        os.chdir(self.whereIStart)

    def runCommands(self,commands):
        """run commands in OpenFOAMCase dir"""
        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        # -- runCommand
        print("I will execute commands in %s:"%self.dir)
        for commandI in range(len(commands)):
            command = commands[commandI]
            print("\t%d. %s"%(commandI+1,command))
            os.system(command)
        
        # -- move back where I start
        os.chdir(self.whereIStart)
    
    def updateTimes(self):
        """update the OpenFOAMCase.times variable"""

        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        lstDir = os.listdir()
        flLstDir = [float(flDir) for flDir in lstDir if isFloat(flDir)]
        self.times = flLstDir
        self.latestTime = max(flLstDir)
        print("Found OpenFOAMCase times are:\n\t%s,\nlatestTime:\n\t%g"%(str(self.times),self.latestTime))

        # -- move back where I start
        os.chdir(self.whereIStart)  
        
    def updateTimesParallel(self):
        """update the OpenFOAMCase.times variable"""

        # -- move to OpenFOAMCase directory 
        os.chdir(self.dir)

        lstDir = os.listdir('processor0')
        flLstDir = [float(flDir) for flDir in lstDir if isFloat(flDir)]
        self.parTimes = flLstDir
        self.latestParTime = max(flLstDir)
        print("Found OpenFOAMCase parallel times are:\n\t%s,\nlatestTime:\n\t%g"%(str(self.parTimes),self.latestParTime))

        # -- move back where I start
        os.chdir(self.whereIStart)  
        
    
