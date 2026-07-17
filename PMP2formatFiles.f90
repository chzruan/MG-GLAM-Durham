Format of data for Bolshoi and GLAM simulations. 
  
File Names: each snapshot consists of a small header file and a number of large data files.
           The header file has information on a given snapshot and on the run. The header
           gives the step number, the expansion parameter and so on. The data
           files provide information on the coordinates and velocities
           of particles. The time-step of the snapshot is coded into
           the file names by using 4 digits of the time-step number.
           For example, the time-step number 17 gives string 0017, which
           is incorporated into the file names. The header filename
           starts with 'PMcrd.' followed by the four digits of the
           time-step followed by ending '.DAT'. The same convention is
           used for the data files. Their names start with 'PMcrs'
           followed by the file number (0, 1, 2 ...). After that the 4 digits of
           the step-number and the '.DAT' ending. For example, file names for
           the time-step 17 may look as:
           PMcrd.0017.DAT, PMcrs0.0017.DAT, PMcrs1.0017.DAT
           The actual number of data files PMcrs depends on the total number of particles.

GLAM-style format of data files.
           


PM-style format of data files.
