//  RealSolutionType.cpp
//
//  Author:
//       Juan J. Durillo <durillo@lcc.uma.es>
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
#include "gobel.h"
#include "RealSolutionType.h"
#include <cstddef>


/**
 * Constructor
 * @param problem
 */
RealSolutionType::RealSolutionType(Problem *problem)
: SolutionType(problem) { }


/**
 * Creates the variables of the solution
 * @param decisionVariables
 */
Variable **RealSolutionType::createVariables() {
  int i;
  double sum = 0,sum1;
  double low_T = 0;
  double tem;
  Variable **variables = new Variable*[problem_->getNumberOfVariables()]; //malloc(sizeof(Real) * problem->getNumberOfVariables());
  if (variables ==  NULL) {
    cout << "Error grave: Impossible to reserve memory for variable type" << endl;
    exit(-1);
  }

  for (i = 0; i < problem_->getNumberOfVariables(); i++) {
    variables[i] = new Real(problem_->getLowerLimit(i),problem_->getUpperLimit(i));
    sum += variables[i]->getValue();
    
    low_T += variables[i]->getLowerBound();
  }
  if (sum > T_MAX_RESOURCE)
  {
      sum1 = 0;
      double d= PseudoRandom::randDouble(0, 1);
      for (i = 0; i < problem_->getNumberOfVariables(); i++)
      {    
      
          tem= variables[i]->getLowerBound()+((variables[i]->getValue()-variables[i]->getLowerBound())*(T_MAX_RESOURCE - low_T) * d / (sum - low_T));
          sum1 += tem;
          variables[i]->setValue(tem);
      }
      
      if (sum1 > T_MAX_RESOURCE)
          cout << "error222" << endl;

  }
  return variables;
} // createVariables
