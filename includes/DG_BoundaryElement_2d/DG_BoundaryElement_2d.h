/**
 * @file DG_BoundaryElement_2d.h
 * @Synopsis  This is the header file for the class which defines the DG_BoundaryElement_2d class. This class  inherits from DG_Element_2d 
 * and adds additional support required for Boundary Elements.
 * @author Krishna Dutt
 * @version 1.0
 * @date 2017-06-19
 */



/* Copyright (c) 2016, Shivasubramanian Gopalakrishnan , Kaushik Kulkarni and Krishna Dutt.
  All rights reserved.

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

     -Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
     -Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS  **AS IS** AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     */

#ifndef DG_BOUNDARYELEMENT_2D_H_
#define DG_BOUNDARYELEMENT_2D_H_

#include "../Utilities/HeaderFiles.h"
#include "../DG_Element_2d/DG_Element_2d.h"

using namespace std;

class DG_BoundaryElement_2d : public DG_Element_2d {
protected:

        map<string, function<void(void)> > TopBoundary;
        map<string, function<void(void)> > BottomBoundary;
        map<string, function<void(void)> > LeftBoundary;
        map<string, function<void(void)> > RightBoundary;

public:

        DG_BoundaryElement_2d(int _N, double x1, double y1, double x2, double y2);
        ~DG_BoundaryElement_2d();

        void assignBoundary(string v, string type, string b);
        void setBoundaryValue(string v, string b);



};

#endif
