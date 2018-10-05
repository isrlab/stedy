[![Documentation](https://img.shields.io/badge/documentation-user_guide-brightgreen.svg)](UserGuide.md)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
# Summary

STEDY is a software package based on MATLAB to enable researchers to simulate dynamics of tensegrity structures. We have developed a Lagrangian formulation for deriving the algebraic differential equations governing the dynamics of classical tensegrity systems. The framework is general enough to allow modeling of general multi-body systems with actuated joints. Furthermore, with the help of the direct correction method that minimizes both geometric and energy constraint violations, errors that arise from numerical integration are corrected on the position and velocity levels, thereby improving simulation accuracy.

The software is aimed at researchers familiar with tensegrity structures, but can be adopted by users from other backgrounds as well. However, users are required to have a preliminary understanding of MATLAB's interface in order to make the best use of this package.

## Dependencies
* [MATLAB Optimization toolbox]

## Startup:
1. Download the STEDY folder from Github.

2. Open MATLAB and make STEDY your working directory.

3. Run setup.m.

4. In case an error message is displayed, follow the instructions below:
    ```MATLAB
    Editing ode solver failed. Please edit ODE45m manually.
    ```
    You will find a copy of *ode45* in your main folder called *ode45m*. Edit it to add the following line exactly in the position shown below. To clarify, in the context of ode45, constraint correction has been implemented only after ensuring that the weighted error (obtained after advancing one step) is no more than the tolerance *rtol*.
    **NOTE**: The user is required to add **ONLY** the line calling the ConstraintCorrection function to ode45m.
    ```matlab
    % Accept the solution only if the weighted error is no more than the
    % tolerance rtol.  Estimate an h that will yield an error of rtol on
    % the next step or the next try at taking this step, as the case may be,
    % and use 0.8 of this value to avoid failures.

    if err > rtol % Failed step
      ...

          if nofailed
            nofailed = false;
            if NNrejectStep
              absh = max(hmin, 0.5*absh);
            else
              absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
            end
          else
            absh = max(hmin, 0.5 * absh);
          end
          h = tdir * absh;
          done = false;

    else  % (Successful step)

        % ADD ONLY THE FOLLOWING LINE CALLING THE FUNCTION
        ynew = ConstraintCorrection(ynew,odeArgs{1},tnew);

          NNreset_f7 = false;
          if nonNegative && any(ynew(idxNonNegative)<0)

      ...
    end
    ```
5. After editing, rename the function at the top to ode45m to avoid a warning.
    ```matlab
    function varargout = ode45m(ode,tspan,y0,options,varargin)
    ```

6. To avoid having to run setup.m everytime you run MATLAB, go to *Set Path* in the Home tab and add the corresponding location of *stedy/Main* to MATLAB's path.

Follow *stedy*'s [User Guide] to get started.
### References
1. Shao-Chen Hsu, Vaishnav Tadiparthi and Raktim Bhattacharya, "A Lagrangian Formulation for Constrained Multibody Dynamics in Tensegrity Systems", Manuscript submitted for publication.

    Please contact the authors at addyhsu@tamu.edu, vaishnavtv@tamu.edu, or raktim@tamu.edu for a copy of the submitted paper.

[User Guide]: (UserGuide.md)
[MATLAB Optimization toolbox]: https://www.mathworks.com/products/optimization.html
