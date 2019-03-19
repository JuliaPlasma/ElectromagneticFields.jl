var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#ElectromagneticFields.jl-1",
    "page": "Home",
    "title": "ElectromagneticFields.jl",
    "category": "section",
    "text": "Common Interface for Electromagnetic Fields(Image: Build Status) (Image: Coverage Status) (Image: codecov)"
},

{
    "location": "#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "Antoine J. Cerfon, Jeffrey P. Freidberg. \"One size fits all\" analytic solutions to the Grad–Shafranov equation. Physics of Plasmas 17 (3), 032502."
},

{
    "location": "#License-1",
    "page": "Home",
    "title": "License",
    "category": "section",
    "text": "Copyright (c) 2017- Michael Kraus <michael.kraus@ipp.mpg.de>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
},

{
    "location": "modules/#",
    "page": "Modules",
    "title": "Modules",
    "category": "page",
    "text": ""
},

{
    "location": "modules/#ElectromagneticFields.ABC",
    "page": "Modules",
    "title": "ElectromagneticFields.ABC",
    "category": "type",
    "text": "ABC equilibrium in (x,y,z) coordinates.\n\nParameters:     A:     B:     C:\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.AxisymmetricTokamakCartesian",
    "page": "Modules",
    "title": "ElectromagneticFields.AxisymmetricTokamakCartesian",
    "category": "type",
    "text": "Axisymmetric tokamak equilibrium in (x,y,z) coordinates.\n\nParameters:     R₀: position of magnetic axis     B₀: B-field at magnetic axis     q:  safety factor\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.AxisymmetricTokamakCylindrical",
    "page": "Modules",
    "title": "ElectromagneticFields.AxisymmetricTokamakCylindrical",
    "category": "type",
    "text": "Axisymmetric tokamak equilibrium in (R,Z,ϕ) coordinates.\n\nParameters:     R₀: position of magnetic axis     B₀: B-field at magnetic axis     q:  safety factor\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.AxisymmetricTokamakToroidal",
    "page": "Modules",
    "title": "ElectromagneticFields.AxisymmetricTokamakToroidal",
    "category": "type",
    "text": "Axisymmetric tokamak equilibrium in (r,θ,ϕ) coordinates.\n\nParameters:     R₀: position of magnetic axis     B₀: B-field at magnetic axis     q:  safety factor\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.Solovev",
    "page": "Modules",
    "title": "ElectromagneticFields.Solovev",
    "category": "type",
    "text": "Axisymmetric Solov\'ev equilibra in (R/R₀,Z/R₀,phi) coordinates. Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010.\n\nParameters:     R₀: position of magnetic axis     B₀: B-field at magnetic axis     ϵ:  inverse aspect ratio     κ:  elongation     δ:  triangularity     a:  free constant, determined to match a given beta value\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.SolovevQuadratic",
    "page": "Modules",
    "title": "ElectromagneticFields.SolovevQuadratic",
    "category": "type",
    "text": "Quadratic Solov\'ev equilibrium in (R,Z,phi) coordinates. Based on McCarthy, Physics of Plasmas 6, 3554, 1999.\n\nParameters:     R₀: position of magnetic axis     B₀: B-field at magnetic axis     a,b: free constants\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.SolovevXpoint",
    "page": "Modules",
    "title": "ElectromagneticFields.SolovevXpoint",
    "category": "type",
    "text": "Axisymmetric Solov\'ev equilibra with X-point in (R/R₀,Z/R₀,phi) coordinates. Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010.\n\nParameters:     R₀: position of magnetic axis     B₀: B-field at magnetic axis     ϵ:  inverse aspect ratio     κ:  elongation     δ:  triangularity     a:  free constant, determined to match a given beta value     xₛₑₚ: x position of the X point     yₛₑₚ: y position of the X point\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.SymmetricQuadratic",
    "page": "Modules",
    "title": "ElectromagneticFields.SymmetricQuadratic",
    "category": "type",
    "text": "Symmetric quadratic equilibrium in (x,y,z) coordinates.\n\nB(xyz) = B_0  beginpmatrix\n0 \n0 \n1 + x^2 + y^2 \nendpmatrix\n\nParameters: B₀\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.ThetaPinch",
    "page": "Modules",
    "title": "ElectromagneticFields.ThetaPinch",
    "category": "type",
    "text": "θ-pinch equilibrium in (x,y,z) coordinates.\n\nParameters:     B₀: B-field at magnetic axis\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.load_equilibrium",
    "page": "Modules",
    "title": "ElectromagneticFields.load_equilibrium",
    "category": "function",
    "text": "Evaluate functions for evaluating analytic equilibria.\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.generate_equilibrium_code",
    "page": "Modules",
    "title": "ElectromagneticFields.generate_equilibrium_code",
    "category": "function",
    "text": "Generate code for evaluating analytic equilibria.\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.generate_equilibrium_functions-Tuple{AnalyticEquilibrium,AnalyticPerturbation}",
    "page": "Modules",
    "title": "ElectromagneticFields.generate_equilibrium_functions",
    "category": "method",
    "text": "Generate functions for evaluating analytic equilibria.\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.get_vector_component-Tuple{Any,Any,Any}",
    "page": "Modules",
    "title": "ElectromagneticFields.get_vector_component",
    "category": "method",
    "text": "Returns the i-th component of the vector corresponding to the one-form α\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.hodge²¹-NTuple{4,Any}",
    "page": "Modules",
    "title": "ElectromagneticFields.hodge²¹",
    "category": "method",
    "text": "Returns the m-th component of the one-form corresponding to the two-form β\n\n\n\n\n\n"
},

{
    "location": "modules/#ElectromagneticFields.jl-–-Modules-1",
    "page": "Modules",
    "title": "ElectromagneticFields.jl – Modules",
    "category": "section",
    "text": "Modules = [ElectromagneticFields]\nOrder   = [:constant, :type, :macro, :function]"
},

]}
