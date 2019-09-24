Regarding parameter files:

1. Currently folder eROSITA contains rmf and arf. Another folder can be made for a different instrument. 

2. File cosmoparams.dat consists of parameters [H0	Ob0	Om0	Ode0	Tcmb0	m_mu] == [Hubble constant, Omega baryon (z=0),
   Omega matter (z=0), Omega DE (z=0), Temperature of CMB (z=0), total mass of neutrinos]. NOTE: Only for constant DE parameter
   w = -1.
3. File mfparams.dat consists of mass function parameters [hms	delta_h	delsel	Mmin	Mmax	dlog10m] == [Selection of mass 
   function, delta_h, selection of delta w.r.t, Minimun mass in sample, Maximum mass in sample, dlog10m (log space bin size)]. 
   NOTE: Use hms = [0, 1, 2, 3, 4] == ['SMT', 'Jenkins', 'Warren', 'Tinker08', 'Tinker10'] & delsel = [0, 1] == ['mean', 'crit']

4. File tfparams.dat consists of transfer function parameters [n_s	sigma_8	lnk_min	lnk_max	dlnk	trans_model_sel] == 
   [primordial spectral index of scalar fluctuations, amplitude of the (linear) power spectrum on the scale of 8 h-1 Mpc, 
   log(Minimum wave vector k), log(Maximum wave vector k), log(bin size of wave vector k)]
   NOTE: Use trans_model_sel = [0, 1, 2, 3] == ['BBKS', 'BondEfs', 'CAMB', 'EH']
   
5. File sampleparams.dat consists of [z_min	z_max	n_z	type] == [Minimum redshift of sample, Maximum redshift of sample, 
   number of redshift bins, type of sample]. NOTE: Use Type = [0, 1] == [linear, logspace]
   
6. File sr_file.dat describes the options of scaling relations. NOTE: Use the code defined in the function in next file 
   srparams.dat
   
7. File srparams.dat consists of all information regarding the scaling relations. Row sr_file	can be used to 
   chose different scaling relations (choice code in sr_file.dat). Row sig_0 is to enter sigma_10 of the different properties 
   at z = 0. Row z_ev_sig_0	can be used to set redshift & mass evolution of sigma (Use 1 for sigma(M,z) and 0 for const sigma).
   Rows norm, slope, z_evo are for the wanted scaling relations Normalization, slope and redshift evolution. The matrix cor 
   describes the correlation coefficient matrix between the different properties. 

8. File xspecparams.dat consists of [abund	El	Eh	Z	ni	expt] == [For setting solar abundance table, min(energy band), 
  max(energy band), solar abundance, reduce original size of table by ni for creating interpolation table, exposure time]
