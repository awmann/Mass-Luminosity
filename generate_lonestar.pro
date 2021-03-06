;; this will generate a set of lonestar scripts

PRO generate_lonestar

  longleaf = 1
  smallstep =   10000
  largestep =  750000
  ;;threads = 72

  close,/all
  openw,13,'batch.sh'
  for err = 1,1 do begin
     for fehon = 0,0 do begin
        for order = 3,3 do begin
           if err ge 1 then erradd = 1 else erradd = 0
           params = order+fehon+erradd+1
           a = strtrim(string(params,format="(I2)"),2)
           b = string(fehon,format="(I1)")
           c = string(err,format="(I1)")
           openw,12,'ML_batch_'+a+'_'+b+'_'+c+'.txt'
           printf,12,'#!/bin/bash'
           printf,12,'#SBATCH -J ML'+a+'_'+b+'_'+c+'           # job name'
           printf,12,'#SBATCH -o ML'+a+'_'+b+'_'+c+'.o%j	# Name of the output file (eg. myMPI.oJobID)'
           printf,12,'#SBATCH -N  4	# number of nodes requested'
           printf,12,'#SBATCH -n  72	# total number of mpi tasks requested'
           if longleaf eq 1 then printf,12,'#SBATCH --mem=10g'
           if longleaf eq 0 then printf,12,'#SBATCH -p normal	# queue (partition) -- normal, development, etc.'
           if longleaf eq 1 then printf,12,'#SBATCH -t 1-12:00:00' else printf,12,'#SBATCH -t 36:00:00	# run time (hh:mm:ss)'
           printf,12,'#SBATCH --mail-user=mann.andrew.w@gmail.com'
           printf,12,'#SBATCH --mail-type=begin   # email me when the job starts'
           printf,12,'#SBATCH --mail-type=end     # email me when the job finishes'
           if longleaf eq 1 then printf,12,'module load python/2.7.12' else printf,12,'module load python/2.7.11'
           printf,12,'python MK-mass_lonestar.py '+a+' '+b+' '+c+' '+strtrim(string(smallstep,format="(I7)"),2)+' '+strtrim(string(largestep,format="(I8)"),2)+' 72'
           close,12
           printf,13,'sbatch ML_batch_'+a+'_'+b+'_'+c+'.txt'
        endfor
     endfor
  endfor
  close,/all

END 
