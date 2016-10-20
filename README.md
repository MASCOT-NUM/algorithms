# template
This repository contains several empty template files for many languages (see http://www.gdr-mascotnum.fr/template.html).

Basics:

- The template is provided for many languages, so I/O objects have to remain in common types (numeric arrays/matrices, strings).
- Four methods/functions are required:
   - 'MyAlgorithm', which allows to initialize the environment, parse options, load dependencies
   - 'getInitialDesign', which returns the first design of experiments
   - 'getNextDesign', which returns the new experiments to perform, or null/empty object if algorithm ends.
   - 'displayResults', which returns a textual/html analysis, possibly including images (in base64 objects)
- Alongside with these requirements, many functions/methods may be appended
- Extended information may be provided in header commented lines (information, algorithm options, dependencies, authors, ...)


These implementations are designed to follow such a workflow:

1. setup the algorithm options according to 'option:' header line
2. algorithm = MyAlgorithm(options) to return the algorithm with its options initialized
3. X0 = getInitialDesign(algorithm,d) to return a first list of experiments, as a matrix of 'd' columns (say 'X0')
4. (external) perform simulations to get each response (say 'Y0') for each matrix X0 line
5. 'X1 = getNextDesign(algorithm,X0,Y0)' to return a new list of experiments, as a new matrix of 'd' columns (say 'X1')
6. (external) perform simulations to get each response (say 'Y1') for each matrix X1 line
7. [optional] 'displayResults(algorithm,X1,Y1)' to return a temporary analysis of results, as an HTML string
8. ...
9. 'Xj = getNextDesign(algorithm,Xi,Yi)' to return a new list of experiments, as a new matrix of 'd' columns (say 'Xj')
10. (external) perform simulations to get each response (say 'Yj') for each matrix 'Xj' line
11. ...
12. displayResults(algorithm,Xn,Yn) to return a final analysis of results, as an HTML string (possibly including base64 files)


![Analytics](https://ga-beacon.appspot.com/UA-109580-20/MASCOT-NUM/template)