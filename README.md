Automated Computation of Spectral Submanifolds 
for Nonlinear Modal Analysis
-----------------------------------------------------------------------------  
Sten Ponsioen & George Haller (ETH Zurich) 

License
-----------------------------------------------------------------------------  
This software is made public for research use only. It may be modified and redistributed
under the terms of the GNU General Public License. 

SSMtool
-----------------------------------------------------------------------------  
SSMtool is an automated computational algorithm for computing two-dimensional spectral submanifolds (SSMs) [2]. 
The algorithm can handle non-conservative mechanical systems of arbitrary (finite) degrees of freedom. We used 
the parameterization method, allowing us to construct the SSMs, their reduced dynamics and corresponding backbone 
curves up to any required order of precision. 

Addendum: Isolas (17 December 2018)
-----------------------------------------------------------------------------  
We have improved the core of SSMtool to handle systems with time-periodic forcing. We have used the exact reduced dynamics 
on two-dimensional time-periodic spectral submanifolds (SSMs) to extract forced-response curves (FRCs) and predict 
isolas in arbitrary multi-degree-of-freedom mechanical systems without performing costly numerical simulations [3].
We are currently implementing this improvement into SSMtool 2.0, which will be released in 2019. 




Citation
-----------------------------------------------------------------------------  
Please cite [1] if you use SSMtool in your own work.

References
-----------------------------------------------------------------------------  
[1] S. Ponsioen, T. Pedergnana & G. Haller, "Automated computation of 
    autonomous spectral submanifolds for nonlinear modal analysis",
    J. Sound Vib. 420 (2018) 269-295. 
	
[2] G. Haller, and S. Ponsioen, "Nonlinear normal modes and spectral 
    submanifolds: existence, uniqueness and use in model reduction",
    Nonlinear Dynamics, 86(3):1493–1534, (2016).

[3] S. Ponsioen, T. Pedergnana & G. Haller, "Analytic Prediction of 
    Isolated Forced Response Curves from Spectral Submanifolds",
    submitted, (2018).

Installation notes 
-----------------------------------------------------------------------------
Tested on MATLAB R2016b, R2017b and R2018b.

1) After you unzipped the files to mydir, 
   put the Current Directory in MATLAB to mydir. 

2) In the MATLAB command prompt,
   type “SSM” to open the graphical user interface.


NOTE: This code may be improved and subject to several changes. Therefore, we suggest to visit this 
page and check if the version you downloaded is up to date.  


Maintained by Sten Ponsioen,
stenp at ethz dot ch
February 26, 2019.
