//  ProblemFactory.cpp
//
//  Author:
//       Esteban López-Camacho <esteban@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "stdafx.h"
#include "ProblemFactory.h"


/**
 * Problem factory
 * @param name : Name of the problem
 */
Problem * ProblemFactory::getProblem(char * name) {
  return getProblem(name, 0, NULL);
}


/**
 * Problem factory which uses the same arguments as the main from where
 * is called (minimum two arguments)
 * @param argc : Number of arguments
 * @param argv : Array of arguments
 */
Problem * ProblemFactory::getProblem(int argc, char ** argv) {
  if (argc==2) {
    return getProblem(argv[1], 0, NULL);
  } else if (argc>2) {
    char * argv1 = argv[1];
    for (int i=0; i<argc-2; i++) {
      argv[i] = argv[i+2];
    }
    return getProblem(argv1, argc-2, argv);
  } else {
    cerr << "Too few arguments to build a problem.";
    exit(-1);
  }
}


/**
 * Problem factory with some optional parameters to be used to construct the
 * problem
 * @param name : Name of the problem
 * @param argc : Number of parameters
 * @param argv : Array of parameters
 */
Problem * ProblemFactory::getProblem(char * name, int argc, char ** argv) {

   if (strcmp(name, "Tanaka")==0) { // Tanaka
  
      return new Tanaka("Real");
    
  }

}

