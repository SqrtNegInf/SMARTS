c
c screen.f -- screen compound to identify scaffolds needing correction
c NB: Some screens are just 'warnings', and add no ClogP correction value
c
c Scaffold screens consist of a mixture of custom code and SMARTS targets. The trade-off
c is that while SMARTS targets are relatively easy to write, they are slow at runtime. Custom
c code is generally 10x or more faster, but harder to write/maintain.
c
c To reduce the runtime penalty for SMARTS, a quick analysis is done for the input compound,
c tallying atoms, bonds, rings (sizes and counts), basic structural features, etc.  This is 
c used to avoid running a lengthy SMARTS search when it cannot possibly succeed (say the target
c contains two nitrogens, but the compound has only one).
c
c A good comparison between the two approaches can be seen in the current SMARTS-based 
c 'phenothiazine' code (around line 500) and the previously-used custom code, left commented 
c just after.
c  
      subroutine  screen( warn, iflag )
      implicit none
      include '$BB_ROOT/h/maxima.inc'
      include '$BB_ROOT/h/clogp.inc'
      include '$BB_ROOT/h/hbond.inc'
      include '$BB_ROOT/h/genie.inc'
c
      real      dummy, dsink, total, value, hbond_factor, hbondval, me_ether, cy_bz_ester, clox_flox,
     *          bz_tunnel, phen_val, ph_thiaz1, ph_thiaz2, pyr_dpalky, pyr_dparom, amide_nh,
     *          pyr_dpnone, pyr_dpnfix, pl_amide, tri_amide, lactone_dibenz, purine_n_chain, 
     *          purine_n_ring, morph, aziridinyl, quinone, azaindole, arom_ccf3,
     *          styrene_mono, styrene_di, styrene_cinn, piper_phenyl, glucopyranoside, 
     *          oxazol_peptide, extend_arom, pyr_ortho, gambogic_acid, peptide, piper_para,
     *          crown_ether, crown_benzo, clogp_max, clogp_min, ccap
      integer   cyal, cyar, cymax, cymin, cy3, cy5, cy6, cyar5, cyal5, cyar6, cyal6, fatom, 
     *          al_s, al_n, al_o, al_or, al_nr, ar_n, ar_nh, ar_o, al_c, at_bz, at_st, al_carb, 
     *          ar_carb, ar_c, at_f, at_cl, at_br, at_halo, db_cc, db_ccr, ar_cn, al_n3,
     *          al_eo, ar_eo, al_amide, ar_amide, at_cf3, ar_s, ch3, s184, cnn, at_oh, al_nhr
      integer   iflag, i, j, jj, k, kk, l, bzt_fr, bzt_arom, bzt_cy, bzt_n, c1, c2, c3, c4, n1, n2, n3
      character cvar40*40
      logical   warn, got_carbonyl, got_n_h, bzt_no2, seen_cy(mxcy), pyrazole, pyrazine, piperidine, piperazine,
     *          got_n_alkyl, got_n_arom, got_f, pass1, sok, fused_ar, arar
      character envi(3)*10, cclass*10, ctype*6, bzt_f*80
c
      logical   fc_has, debug, debugg, debugv, r186, r187
      integer   list_i(3,2), b, tmp1, tmp2
      data pass1 /.true./
      data debug /.false./
c
c these flags can be enabled to collect information for performance auditing
c
      data debugg /.false./
      data debugv /.false./
c
c rule 185 blocks all scaffold corrections
c rule 186 blocks subset: all non-SMARTS 
c rule 187 blocks subset: all SMARTS
c
      if (.not. audit_status(185)) return
      r186 = audit_status(186)
      r187 = audit_status(187)
c       print *, 'r187', r187
c
c Read values from constant database CNSTDB on 1st pass
c
      if (pass1) then
         call cnstdb( 'FIND', 'ME-ETHER', me_ether )
         call cnstdb( 'FIND', 'CYBZESTER', cy_bz_ester )
         call cnstdb( 'FIND', 'BZ-TUNNEL', bz_tunnel )
         call cnstdb( 'FIND', 'PHTHIAZ1', ph_thiaz1 )
         call cnstdb( 'FIND', 'PHTHIAZ2', ph_thiaz2 )
         call cnstdb( 'FIND', 'PYRDPNONE', pyr_dpnone )
         call cnstdb( 'FIND', 'PYRDPNFIX', pyr_dpnfix )
         call cnstdb( 'FIND', 'PYRDPAROM', pyr_dparom )
         call cnstdb( 'FIND', 'PYRDPALKY', pyr_dpalky )
         call cnstdb( 'FIND', 'PL_AMIDE', pl_amide )
         call cnstdb( 'FIND', 'HBONDVAL', hbondval )
         call cnstdb( 'FIND', 'PURINE-N-RING', purine_n_ring )
         call cnstdb( 'FIND', 'PURINE-N-CHAIN', purine_n_chain )
         call cnstdb( 'FIND', 'TRI-AMIDE', tri_amide )
         call cnstdb( 'FIND', 'LACTONE-DIBENZ', lactone_dibenz )
         call cnstdb( 'FIND', 'AMIDE-NH', amide_nh )
         call cnstdb( 'FIND', 'CLOX-FLOX', clox_flox )
         call cnstdb( 'FIND', 'MORPHOLINE', morph )
         call cnstdb( 'FIND', 'QUINONE-AMINE', quinone )
         call cnstdb( 'FIND', 'AZIRIDINYL', aziridinyl )
         call cnstdb( 'FIND', 'AZAINDOLE', azaindole )
         call cnstdb( 'FIND', 'AROM-CCF3', arom_ccf3 )
         call cnstdb( 'FIND', 'STYRENE-MONO', styrene_mono )
         call cnstdb( 'FIND', 'STYRENE-DI', styrene_di )
         call cnstdb( 'FIND', 'STYRENE-CINN', styrene_cinn )
         call cnstdb( 'FIND', 'PIPER-PHENYL', piper_phenyl )
         call cnstdb( 'FIND', 'GLUCO-PYRANO', glucopyranoside )
         call cnstdb( 'FIND', 'OXAZOL-PEPTIDE', oxazol_peptide )
         call cnstdb( 'FIND', 'EXTEND-AROM', extend_arom )
c        call cnstdb( 'FIND', 'PYRAZINE-CYANO', pyr_cyano )
c        call cnstdb( 'FIND', 'PYRAZINE-HALO', pyr_halo )
c        call cnstdb( 'FIND', 'PYRAZINE-CARBOXY', pyr_carboxy )
         call cnstdb( 'FIND', 'PYRAZINE-ORTHO', pyr_ortho )
         call cnstdb( 'FIND', 'GAMBOGIC-ACID', gambogic_acid )
         call cnstdb( 'FIND', 'PEPTIDE', peptide )
         call cnstdb( 'FIND', 'PIPERIDINE-PARA', piper_para )
         call cnstdb( 'FIND', 'CROWN-ETHER', crown_ether )
         call cnstdb( 'FIND', 'CROWN-BENZO', crown_benzo )
         pass1 = .false.
      endif
c
c contribution class description
c
      cclass = 'Screen'
      ctype  = ' '
c
      iflag = 0
c
c characterize compound in detail, so SMARTS searches can be skipped, for efficiency
c
c ring analysis
c
      arar  = .false.
      cyal  = 0
      cyar  = 0
      cy3   = 0
      cy5   = 0
      cyal5 = 0
      cyar5 = 0
      cy6   = 0
      cyal6 = 0
      cyar6 = 0
      cymax = 0
      cymin = 99 
      pyrazole   = .false.
      pyrazine   = .false.
      piperidine = .false.
      piperazine = .false.
      do 10 i=1,ncycle
         if (cylen(i).gt.cymax) cymax = cylen(i)
         if (cylen(i).lt.cymin) cymin = cylen(i)
         if (cylen(i).eq. 3) cy3   = cy3 + 1
         if (cylen(i).eq. 5) cy5   = cy5 + 1
         if (cylen(i).eq. 5 .and.       cyarom(i)) cyar5   = cyar5 + 1
         if (cylen(i).eq. 5 .and. .not. cyarom(i)) cyal5   = cyal5 + 1
         if (cylen(i).eq. 6 .and.       cyarom(i)) cyar6   = cyar6 + 1
         if (cylen(i).eq. 6 .and. .not. cyarom(i)) cyal6   = cyal6 + 1
         if (cylen(i).eq. 6) cy6   = cy6 + 1
         if (.not.cyarom(i)) cyal  = cyal + 1
         if (     cyarom(i)) cyar  = cyar + 1
         if (cyarom(i) .and. cylen(i).eq.5) then
            n1 = 0
            n2 = 0
            n3 = 0
            do 7 j=1,cylen(i)
               if (atsymb(cycle(i,j)).eq.'n ') then
                  n3 = n3 + 1
                  if (n1.eq.0) n1 = j
                  if (n1.ne.0) n2 = j
               endif
  7         continue
            if (n3.eq.2 .and. (n2-n1).eq.1 .or. (n1.eq.1 .and. n2.eq.5)) pyrazole = .true.
         elseif (cyarom(i) .and. cylen(i).eq.6) then
            c1 = 0
            c2 = 0
            n1 = 0
            n2 = 0
            do 8 j=1,cylen(i)
               if (atsymb(cycle(i,j)).eq.'n ') then
                  c1 = c1 + 1
                  if (c1.eq.1) n1 = j
                  if (c1.eq.2) n2 = j
               endif
               if (atsymb(cycle(i,j)).eq.'c ') c2 = c2 + 1
  8         continue
            if (.not. (c1.eq.2 .and. c2.eq.4)) goto 10
            if (c1.eq.2 .and. c2.eq.4 .and. ( n2-n1.eq.-3 .or. n2-n1.eq.3 ) ) pyrazine = .true.
         elseif (.not.cyarom(i) .and. cylen(i).eq.6) then
            if (atsymb(cycle(i,1)).eq.'N ' .and. atsymb(cycle(i,4)).eq.'N ') piperazine = .true.
            if (atsymb(cycle(i,2)).eq.'N ' .and. atsymb(cycle(i,5)).eq.'N ') piperazine = .true.
            if (atsymb(cycle(i,3)).eq.'N ' .and. atsymb(cycle(i,6)).eq.'N ') piperazine = .true.
            n1 = 0
            c1 = 0
            do 9 j=1,cylen(i)
               if (atsymb(cycle(i,j)).eq.'C ') c1 = c1 + 1
               if (atsymb(cycle(i,j)).eq.'N ') n1 = n1 + 1
  9         continue
            if (c1.eq.5 .and. n1.eq.1) piperidine = .true.
         endif
 10   continue
c
c atom analysis
c
      fatom = 0
      ch3   = 0
      cnn   = 0
      al_s  = 0
      s184  = 0
      ar_s  = 0
      al_n  = 0
      al_n3 = 0
      ar_n  = 0
      ar_nh = 0
      al_o  = 0
      al_or = 0
      al_nr = 0
      al_nhr= 0
      ar_o  = 0
      ar_cn = 0
      al_c  = 0
      ar_c  = 0
      al_eo = 0
      ar_eo = 0
      at_f  = 0
      at_oh = 0
      at_cl = 0
      at_br = 0
      at_bz = 0
      at_st = 0
      al_carb = 0
      ar_carb = 0
      at_cf3  = 0
      al_amide   = 0
      ar_amide   = 0
      db_cc   = 0
      db_ccr  = 0
      fused_ar = .false.
      do 100 i=1,n
         if (  cmap(i).gt.1   ) fatom = fatom + 1
         if (atsymb(i).eq.'S ') al_s  = al_s  + 1
         if (atsymb(i).eq.'s ') ar_s  = ar_s  + 1
         if (atsymb(i).eq.'n' ) ar_n  = ar_n  + 1
         if (atsymb(i).eq.'N ') al_n  = al_n  + 1
         if (atsymb(i).eq.'o' ) ar_o  = ar_o  + 1
         if (atsymb(i).eq.'O ') al_o  = al_o  + 1
         if (atsymb(i).eq.'C ') al_c  = al_c  + 1
         if (atsymb(i).eq.'C ' .and. hcount(i).eq.3) ch3 = ch3 + 1
         if (atsymb(i).eq.'Cl') at_cl = at_cl + 1
         if (atsymb(i).eq.'Br') at_br = at_br + 1
         if (atsymb(i).eq.'F ') at_f  = at_f  + 1
         if (isotyp(imap(i)).eq.'Z') at_bz = at_bz + 1
         if (isotyp(imap(i)).eq.'Y') at_st = at_st + 1
         if (atsymb(i).eq.'n'  .and. hcount(i).eq.1) ar_nh = ar_nh + 1
         if (atsymb(i).eq.'N ' .and. hcount(i).eq.0) al_n3 = al_n3 + 1
         if (atsymb(i).eq.'O ' .and. hcount(i).eq.1) at_oh = at_oh + 1
         if (atsymb(i).eq.'O ' .and. cmap(i).gt.0) al_or  = al_or + 1
         if (atsymb(i).eq.'N ' .and. cmap(i).gt.0) al_nr  = al_nr + 1
         if (atsymb(i).eq.'N ' .and. cmap(i).gt.0 .and. hcount(i).eq.1) al_nhr  = al_nhr + 1
         if (atsymb(i).eq.'S ' .and. fsmile(fmap(i))(1:11).eq.'AS(Z)(=O)=O') s184 = s184 + 1
         if (atsymb(i).eq.'O ' .and. hcount(i).eq.0 .and. ncon(i).eq.2 .and.
     *       atsymb(con(i,1)).eq.'C ' .and. atsymb(con(i,2)).eq.'C ') al_eo = al_eo + 1
         if (atsymb(i).eq.'O ' .and. hcount(i).eq.0 .and. ncon(i).eq.2 .and.
     *       ( (atsymb(con(i,1)).eq.'C ' .and. atsymb(con(i,2)).eq.'c ') .or.
     *         (atsymb(con(i,1)).eq.'c ' .and. atsymb(con(i,2)).eq.'C ') ) ) ar_eo = ar_eo + 1
c
         if (atsymb(i).eq.'C ') then
            do 30 jj=1,ncon(i)
               j = con(i,jj)
               if (atsymb(j).eq.'O ' .and. cmap(j).eq.0 .and. bond(i,jj).eq.2) al_carb = al_carb + 1
               if (j.lt.i) goto 30
               if (atsymb(j).eq.'C ' .and. cmap(j).eq.0 .and. bond(i,jj).eq.2) db_cc  = db_cc  + 1
               if (atsymb(j).eq.'C ' .and.                    bond(i,jj).eq.2) db_ccr = db_ccr + 1
 30         continue
         endif
         if (atsymb(i).eq.'C ' .and. ncon(i).eq.3) then
c        if (atsymb(i).eq.'C ' .and. cmap(i).eq.0 .and. ncon(i).eq.3) then
            c1 = 0
            c2 = 0
            c3 = 0
            c4 = 0
            do 31 jj=1,3
               j = con(i,jj)
               if (atsymb(j).eq.'O ' .and. bond(i,jj).eq.2) c1 = 1
               if (atsymb(j).eq.'N ') c2 = c2 + 1
               if (atsymb(j).eq.'C ') c3 = c3 + 1
               if (atsymb(j).eq.'c ') c4 = c4 + 1
 31         continue
            if (c1.eq.1 .and. c2.eq.1 .and. c3.eq.1) al_amide = al_amide + 1
            if (c1.eq.1 .and. c2.eq.1 .and. c4.eq.1) ar_amide = ar_amide + 1
         endif
         if (atsymb(i).eq.'C ' .and. cmap(i).eq.0 .and. ncon(i).eq.4) then
            c1 = 0
            c2 = 0
            do 32 jj=1,4
               j = con(i,jj)
               if (atsymb(j).eq.'C ') c1 = c1 + 1
               if (atsymb(j).eq.'F ') c2 = c2 + 1
 32         continue
            if (c1.eq.1 .and. c2.eq.3) at_cf3 = at_cf3 + 1
         endif
         if (atsymb(i).eq.'c ') then
            ar_c = ar_c + 1
            c1 = 0
            n1 = 0
            do 40 jj=1,ncon(i)
               j = con(i,jj)
               if (atsymb(j).eq.'n ') n1 = n1 + 1
               if (atsymb(j).eq.'c ' .or. atsymb(j).eq.'n ' .or. atsymb(j).eq.'o ') c1 = c1 + 1
               if (atsymb(j).eq.'O ' .and. cmap(j).eq.0 .and. bond(i,jj).eq.2) ar_carb = ar_carb + 1
               if (atsymb(j).eq.'N '                                         ) ar_cn   = ar_cn   + 1
               if ((atsymb(j).eq.'c ' .or. atsymb(j).eq.'n ') .and. bond(i,jj).eq.1) arar = .true.
 40         continue
            if (n1.ge.2) cnn = n1
            if (c1.ge.3) fused_ar = .true.
         endif
c
 100  continue
      at_halo = at_f + at_cl + at_br
c
c=================================================================================
c corrections start here
c=================================================================================
c
c  steroid and cyclic sugar
c
      if (r186 .and. audit_status(1) .and. n.ge.18 .and. fatom.ge.6 .and. ncycle.ge.4 .and. cyal.ge.3) then
         if (debugg) print *, 'screen: 1 steroids/glycosides (no SMARTS)'
         if (audit_status(1)) call stare(warn)
      endif
c
c planar amide
c
      if (r186 .and. audit_status(10) .and. fatom.ge.4 .and. cyal.ge.2 .and. cyar.ge.1 .and. al_carb.ge.1 .and. cy5.ge.1 .and. cy6.ge.2) then
         if (debugg) print *, 'screen: 10 planar amide (no SMARTS)'
         do 600 i=1,n
            if (atnumb(i).eq.7 .and. 
     *          cmap(i).eq.2 .and.
     *          fsmile(fmap(i))(1:11).eq.'AN(a)C(A)=O'
     *         ) then
               call outlin(cclass,ctype,' Planar amide',' ',pl_amide)
               dummy = total(pl_amide)
            endif
 600     continue
      endif
c
c ketone_biphenyl
c
      if (r186 .and. audit_status(35) .and. cyal.ge.1 .and. cyar.ge.2 .and. fatom.ge.4 .and. al_carb.ge.1 .and. cy6.ge.2) then
         if (debugg) print *, 'screen: 35 ketone_biphenyl (no SMARTS)'
         call ketone_biphenyl(warn)
      endif
c
c cyclopropyl 
c
      if (r186 .and. audit_status(42) .and. cy3.ge.1 .and. nfrag.ge.1) then
         if (debugg) print *, 'screen: 42 cyclopropyl (no SMARTS)'
         call cyclopropyl(warn)
      endif
c
c pyridone dipole factor
c  skip correction for aromatic fused rings
c  added '*nn(*)*' (diazole-n-subst)
c  tweak value based on type of N substitution
c  allow phenyl substitution on N
c  apply correction 6-membered aromatic rings that have both carbonyl fragment & nitrogen
c  discarded Elecias code, completely redone
c
      if (r186 .and. audit_status(4) .and. ar_n.ge.1 .and. ar_carb.ge.1 .and. cy6.ge.1) then
      hbond_factor = 1
      do 710 i=1,ncycle
          if (.not.cyarom(i) .or. 6.ne.cylen(i)) goto 710
          got_carbonyl = .false.
          got_n_h      = .false.
          got_n_alkyl  = .false.
          got_n_arom   = .false.
          got_f        = .false.
          do 709 j=1,6
              k  = cycle(i,j)
              jj = fmap(k)
              if (jj .le. 0) goto 709
              if (fsmile(jj)(1:8) .eq. 'ac(a)=O ') got_carbonyl = .true.
              if (fsmile(jj)(1:7) .eq. 'a[nH]a ' ) got_n_h      = .true.
              if (fsmile(jj)(1:7) .eq. 'An(a)a ' ) got_n_alkyl  = .true.
              if (fsmile(jj)(1:7) .eq. 'Zn(a)a ' ) got_n_alkyl   = .true.  
              if (fsmile(jj)(1:7) .eq. 'an(a)a ' .and. cmap(k).eq.1) got_n_arom = .true.
 709      continue
          if (r187 .and. audit_status(110) .and. got_carbonyl .and. got_n_h) then
              if (debugg) print *, 'screen: 110a pyridone dipole'
              if (debugv) print *, 'vectar: 110 c1[n;H]ccc(=O)c1O'
              vectar(0) = 'c1[n;H]ccc(=O)c1O'
              call search( 0, 'ANYONE', sok )
              value = pyr_dpnone * hbond_factor
              if (nhits .ne. 0) value = value + pyr_dpnfix
              cvar40 = ' Pyridone dipole factor (unsubst)'
              if (hbond_factor.ne.1) cvar40 = ' Pyridone dipole factor (unsubst) [xx%]'
              call outlin(cclass,ctype,cvar40,' ',value)
              dummy = total(value)
          else if (r186 .and. audit_status(111) .and. got_carbonyl .and. got_n_alkyl) then
              if (debugg) print *, 'screen: 110b pyridone dipole (no SMARTS)'
              cvar40 = ' Pyridone dipole factor (alkyl)'
              value = pyr_dpalky * hbond_factor
              if (hbond_factor.ne.1) cvar40 = ' Pyridone dipole factor (alkyl) [xx%]'
              call outlin(cclass,ctype,cvar40,' ',value)
              dummy = total(value)
          else if (r186 .and. audit_status(112) .and. got_carbonyl .and. got_n_arom ) then
              if (debugg) print *, 'screen: 110c pyridone dipole (no SMARTS)'
              cvar40 = ' Pyridone dipole factor (aromatic)'
              value = pyr_dparom * hbond_factor
              if (hbond_factor.ne.1) cvar40 = ' Pyridone dipole factor (arom) [xx%]'
              call outlin(cclass,ctype,cvar40,' ',value)
              dummy = total(value)
          endif
 710  continue
 711  continue
      endif
c
c methoxy ether (aliphatic)
c apply when bonded to CH3 and an aliphatic carbon that is not double-bonded
c
       if (r186 .and. audit_status(14) .and. al_eo.ge.1) then
          do 800 i=1,n
             if (atsymb(i).ne.'O' .or. ncon(i).ne.2) goto 800
             call fc_envv(i,envi,3,list_i,2)
             if (.not. fc_has('2:V-AL',envi,3,b)) goto 800
             tmp1 = con(i,1)
             tmp2 = con(i,2)
             if (1.ne.ncon(con(i,1)) .and. 1.ne.ncon(con(i,2))) goto 800
             if (debugg) print *, 'screen: 14 methyl ether (no SMARTS)'
             call outlin(cclass,ctype,' Methyl ether (aliphatic)',' ',me_ether)
             dummy = total(me_ether)
  800     continue
       endif
c
c cyclic benzyl ester (aromatic env.)
c
      if (r186 .and. audit_status(40) .and. cyar.ge.1 .and. cyal.ge.1 .and. fatom.ge.2 .and. al_carb.ge.1 .and. at_bz.ge.1 .and. al_eo.ge.1) then
        if (debugg) print *, 'screen: 40 cyclic benzyl ester (no SMARTS)'
        do 900 i=1,n
          if (atsymb(i).ne.'C' .or. ncon(i).ne.3 .or. cmap(i).eq.0) goto 900
          if (fmap(i) .le. 0) goto 900
          if (fsmile(fmap(i)) .ne. 'ZOC(a)=O') goto 900
          call outlin(cclass,ctype,' Cyclic benzyl ester',' ',cy_bz_ester)
          dummy = total(cy_bz_ester)
  900   continue
      endif
c
c cyclic benzyl thio-ester (aromatic env.) 
c
      if (r186 .and. audit_status(142) .and. cyar.ge.1 .and. cyal.ge.1 .and. fatom.ge.2 .and. al_carb.ge.1. .and. al_s.ge.1) then
        if (debugg) print *, 'screen: 142 cyclic benzyl thio-ester (no SMARTS)'
        do 901 i=1,n
          if (atsymb(i).ne.'C' .or. ncon(i).ne.3 .or. cmap(i).eq.0) goto 901
          if (fmap(i) .le. 0) goto 901
          if (fsmile(fmap(i)) .ne. 'ZSC(a)=O') goto 901
          call outlin(cclass,ctype,' Cyclic benzyl thio-ester',' ',cy_bz_ester)
          dummy = total(cy_bz_ester)
  901   continue
      endif
c
c benzyl tunneling 
c tunneling through benzyl for constant-valued sigma/rho interaction
c CH2-OH attached to pyridine (or NO2 substituted on same ring)
c
      if (r186 .and. audit_status(41) .and. cyar.ge.1 .and. nfrag.ge.2. .and. al_o.ge.1 .and. at_bz.ge.1) then
        if (debugg) print *, 'screen: 41a benzyl tunneling (no SMARTS)'
        do 910 i=1,ncycle
          seen_cy(i) = .false.
  910   continue
        do 1000 i=1,n
          if (atsymb(i).ne.'C' .or. ncon(i).ne.2 .or. cmap(i).ne.0) goto 1000
          if (fmap(i) .gt. 0) goto 1000
          bzt_fr   = 0
          bzt_arom = 0
          bzt_n    = 0
          bzt_no2  = .false.
          do 990 jj=1,2
             j = con(i,jj)
             if (atsymb(j).eq.'O') then
               bzt_fr  = j
             elseif (atsymb(j).eq.'c') then
               bzt_cy = cymemb(j,1)
               if (cylen(bzt_cy).ne.6) goto 990
                 do 980 kk=1,6
                   k = cycle(bzt_cy,kk)
                   if (atsymb(k).eq.'n') then
                     bzt_n = bzt_n + 1
                   elseif (atsymb(k).eq.'c' .and. ncon(k).eq.3) then
                     do 970 l=1,3
                       if (atsymb(con(k,l)).eq.'N' .and. gsmile(fmap(con(k,l))).eq.'*N(=O)=O') then
                         bzt_no2 = .true.
                       endif
  970                continue
                 endif
  980          continue
               if (bzt_n.eq.1 .or. bzt_no2) bzt_arom = j
            endif
  990     continue
          if (bzt_fr.eq.0 .or. bzt_arom .eq.0) goto 1000
          bzt_f = gsmile(fmap(bzt_fr))
          if (bzt_f(1:2).ne.'*O') goto 1000
          if (0 .ne. index(bzt_f(2:),'*')) goto 1000
          if (seen_cy(bzt_cy)) goto 1000
          seen_cy(bzt_cy) = .true.
          if (bzt_n.eq.1) then
          call outlin(cclass,ctype,' Benzyl tunneling (hetero)',' ',bz_tunnel)
          else
          call outlin(cclass,ctype,' Benzyl tunneling (nitro)',' ',bz_tunnel)
          endif
          dummy = total(bz_tunnel)
 1000   continue
      endif
c
c phenothiazines
c
      if (r187 .and. audit_status(45) .and. cyar.ge.2 .and. cyal.ge.1 .and. fatom.ge.4 .and. al_s.ge.1 .and. al_n3.ge.1) then
         if (debugg) print *, 'screen: 45 phenothiazine'
         if (debugv) print *, 'vectar:  45 [N;H0]2c1ccccc1Sc3ccccc23'
         vectar(0) = '[N;H0]2c1ccccc1Sc3ccccc23'
         call search( 0, 'ANYONE', sok )
         if (nhits .ne. 0) then
            phen_val = ph_thiaz1
            if (0.ne.index(ssmile,'S(=O)')) phen_val = ph_thiaz2
            call outlin(cclass,ctype,' Phenothiazine',' ',phen_val)
            dummy = total(phen_val)
         endif
      endif
c
c XXX phenothiazines XXX obsolete!  
c replaced by SMARTS target code above
c
c     if (audit_status(45)) then
c        if (debugg) print *, 'screen: 45 phenothiazine'
c        phen_s = 0
c        phen_n = 0
c        do 1010 i=1,n
c           if (ncon(i).lt.2 .or. cmap(i).ne.1 .or. fmap(i).le.0) goto 1010
c           if (atsymb(i).eq.'S') then
c              fstr = fsmile(fmap(i))
c              if (fstr.ne.'aSa' .and. fstr.ne.'aS(a)=O' .and. fstr.ne.'aS(a)(=O)=O') goto 1010
c              phen_s = cymemb(i,1)
c              if (fstr.eq.'aSa') then
c                 if (.not. audit_status(xxx)) goto 1010 ! was 108
c                 phen_val = ph_thiaz1
c              elseif (audit_status(xxx)) then ! was 109
c                 phen_val = ph_thiaz2
c              endif
c              goto 1011
c           endif
c1010    continue
c1011    continue
c        if (phen_s .ne. 0) then
c           do 1020 i=1,n
c              if (ncon(i).lt.2 .or. cmap(i).ne.1 .or. fmap(i).le.0) goto 1020
c              if (atsymb(i).eq.'N' .and. hcount(i).eq.0) then
c                 fstr = fsmile(fmap(i))
c                 if (fstr.ne.'aNa' .and. fstr.ne.'AN(a)a' .and. fstr(1:5).ne.'aN(a)') goto 1020
c                 phen_n = cymemb(i,1)
c                 goto 1021
c           endif
c1020    continue
c        endif
c1021    continue
c        if (phen_s.ne.0 .and. phen_s.eq.phen_n) then
c           call outlin(cclass,ctype,' Phenothiazine',' ',phen_val)
c           dummy = total(phen_val)
c        endif
c     endif
c
c SMARTS: tri-amide
c
      if (r187 .and. audit_status(49) .and. al_amide.ge.1 .and. n.ge.18 .and. al_c.ge.6 .and. al_n3.ge.1 .and. al_o.ge.2 .and. cyar.ge.1) then
        if (debugg) print *, 'screen: 49 tri-amide di-ethyl+'
        if (debugv) print *, 'vectar:  49 [A;R0]CCN(CC[A;R0])C([C;R0])=O'
        vectar(0) = '[A;R0]CCN(CC[A;R0])C([C;R0])=O'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Tri-amide (di-ethyl+)',' ',tri_amide)
            dummy = total(tri_amide)
         endif
      endif
c
c SMARTS: purine nucleosides
c separate corrections for -CH-OH on a ring, in a chain
c
      if (r187 .and. audit_status(50) .and. cnn.ge.1 .and. ar_n.ge.3 .and. (ar_c+ar_n).ge.9 .and. fatom.ge.2 .and. 
     *      al_eo.ge.1 .and. al_o.ge.2 .and. cyar6.ge.1 .and. cyar5.ge.1) then
        if (debugg) print *, 'screen: 50a purine nucleoside (ring)'
        if (debugv) print *, 'vectar:  50 [O;H][C,O][C,O][C,O,S][C;R1]n[c;R2]n A'
        vectar(0) = '[O;H][C,O][C,O][C,O,S][C;R1]n[c;R2]n'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Purine nucleoside','Ring',purine_n_ring)
            dummy = total(purine_n_ring)
        elseif (cyal.eq.0) then
            if (debugg) print *, 'screen: 50b purine nucleoside (chain)'
            if (debugv) print *, 'vectar:  50 [O;H][C,O][C,O][C,O,S][C;R0]n[c;R2]n B'
            vectar(0) = '[O;H][C,O][C,O][C,O,S][C;R0]n[c;R2]n'
            call search( 0, 'ANYONE', sok )
            if (nhits .ne. 0) then
                call outlin(cclass,'SMARTS',' Purine nucleoside','Chain',purine_n_chain)
                dummy = total(purine_n_chain)
            endif
        endif
      endif
c
c SMARTS: lactone adjacent to fused ring
c
      if (r187 .and. audit_status(51) .and. at_bz.ge.2 .and. fatom.ge.2 .and. al_eo.ge.1 .and. al_carb.ge.1 .and. cyar6.ge.1 .and. cyal6.ge.1) then
        if (debugg) print *, 'screen: 51 lactone/dibenzyl'
        if (debugv) print *, 'vectar:  51 a[C;R1][C;R1](=O)[O;R1][C;R1]a'
        vectar(0) = 'a[C;R1][C;R1](=O)[O;R1][C;R1]a'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Lactone (dibenzyl)',' ',lactone_dibenz)
            dummy = total(lactone_dibenz)
         endif
      endif
c
c SMARTS: cloxacillin/floxacillin (skip any carboxylic acids)
c
      if (r187 .and. audit_status(53) .and. cyar.ge.1 .and. cy5.ge.1 .and. al_carb.ge.1 .and. ar_n.ge.1 .and. ar_o.ge.1. 
     *     .and. al_o.ge.3 .and. at_bz.ge.1) then
        if (debugg) print *, 'screen: 53a cloxacillin/floxacillin'
        if (debugv) print *, 'vectar:  53 O=C[c;R1][c;R1][C;H2]O A'
        vectar(0) = 'O=C[c;R1][c;R1][C;H2]O'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            if (debugg) print *, 'screen: 53b cloxacillin/floxacillin'
            if (debugv) print *, 'vectar:  53 O=C(O)[c;R1][c;R1][C;H2]O B'
            vectar(0) = 'O=C(O)[c;R1][c;R1][C;H2]O'
c           vectar(0) = 'O=C(O)[c;R1][c;R1][C;H2][O;H]'
            call search( 0, 'ANYONE', sok )
            if (nhits .ne. 0) goto 1100 
            call outlin(cclass,'SMARTS',' Cloxacillin/Floxacillin',' ',clox_flox)
            dummy = total(clox_flox)
         endif
      endif
 1100 continue
c
c SMARTS: Amide - n water bridge
c
      if (r187 .and. audit_status(61) .and. ar_amide.ge.1 .and. cyar.ge.2 .and. fatom.ge.2 .and. ar_n.ge.1 .and. al_carb.ge.1 .and. al_n.ge.1.and. ar_c.ge.7) then
        if (debugg) print *, 'screen: 61b amide - n water bridge (apply correction)'
        if (debugv) print *, 'vectar:  61 [N;H1,H2;R0][C;R0](=O)[c;R1][c;R2][n;R1][c;R1]'
        vectar(0) = '[N;H1,H2;R0][C;R0](=O)[c;R1][c;R2][n;R1][c;R1]'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
           call outlin(cclass,'SMARTS',' Amide/nitrogen water bridge',' ',amide_nh)
           dummy = total(amide_nh)
        endif
      endif
c
c SMARTS: Morpholine
c  include MORPHOLINE itself (and a very few other compounds)
c
      if (r187 .and. audit_status(118) .and. cyal.ge.1 .and. al_nr.ge.1 .and. al_or.ge.1 .and. cy6.ge.1) then
        if (debugg) print *, 'screen: 118a morpholine'
        if (debugv) print *, 'vectar: 118 O1CCNCC1'
        vectar(0) = 'O1CCNCC1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
           if (debugg) print *, 'screen: 118a morpholine'
           if (debugv) print *, 'vectar: 118 O1CCN(C)CC1 A'
           vectar(0) = 'O1CCN(C)CC1'
           call search( 0, 'ANYONE', sok )
           if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Morpholine','N-C',morph)
              dummy = total(morph)
           else if (al_nhr.ge.1) then
              if (debugg) print *, 'screen: 118b morpholine'
              if (debugv) print *, 'vectar: 118 O1CC[N;H]CC1 B'
              vectar(0) = 'O1CC[N;H]CC1'
              call search( 0, 'ANYONE', sok )
              if (nhits .ne. 0) then
                 call outlin(cclass,'SMARTS',' Morpholine','NH',morph)
                 dummy = total(morph)
              endif
           endif
        endif
      endif
c
c SMARTS: benzo/naphtho-quinone with amine ring
c
      if (r187 .and. audit_status(123) .and. al_carb.ge.2 .and. cyal6.ge.1 .and. al_n3.ge.1 .and. 
     *      ( db_ccr.ge.2 .or. (db_ccr.ge.1 .and. cyar.ge.1 .and. fatom.ge.2) ) ) then
        if (debugg) print *, 'screen: 123 benzo/naptho-quinone'
        if (debugv) print *, 'vectar: 123 O=C1[C,c]~[C,c]C(=O)C=C1[N;H0;R1]' 
        vectar(0) = 'O=C1[C,c]~[C,c]C(=O)C=C1[N;H0;R1]'
        call search( 0, 'UNIQUE', sok )
        if (nhits .ne. 0) then
           call outlin(cclass,'SMARTS',' Benzo/naphtho-quinone',' ',nhits*quinone)
           dummy = total(nhits*quinone)
        endif
      endif
c
c SMARTS: Aziridinyl on benzene
c
      if (r187 .and. audit_status(124) .and. cy3.ge.1 .and. cyar6.ge.1 .and. al_n3.ge.1 .and. al_nr.ge.1) then
        if (debugg) print *, 'screen: 124 aziridnyl'
        if (debugv) print *, 'vectar: 124 C1CN1c2ccccc2'
        vectar(0) = 'C1CN1c2ccccc2'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
           call outlin(cclass,'SMARTS',' Aziridinyl',' ',aziridinyl)
           dummy = total(aziridinyl)
        endif
      endif
c
c SMARTS: 7-azaindole
c
      if (r187 .and. audit_status(134) .and. cnn.gt.1 .and. fatom.ge.2 .and. ar_n.ge.2 .and. ar_c.ge.7 .and. cy5.ge.1 .and. cy6.ge.1) then
        if (debugg) print *, 'screen: 134 7-azaindole'
        if (debugv) print *, 'vectar: 134 c1cnc2[n]ccc2c1'
        vectar(0) = 'c1cnc2[n]ccc2c1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
           call outlin(cclass,'SMARTS',' 7-Azaindole',' ',azaindole)
           dummy = total(azaindole)
        endif
      endif
c
c SMARTS: n-C-CF3
c
      if (r187 .and. audit_status(135) .and. at_cf3.ge.1. .and. cyar.ge.1 .and. ar_n.ge.1) then
        if (debugg) print *, 'screen: 135 methyl-CF3'
        if (debugv) print *, 'vectar: 135 n[C;R0]C(F)(F)F'
        vectar(0) = 'n[C;R0]C(F)(F)F'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
           call outlin(cclass,'SMARTS',' Methyl-CF3',' ',arom_ccf3)
           dummy = total(arom_ccf3)
        endif
      endif
c
c SMARTS: 
c
c     if (0.eq.1) then
      if (r187 .and. audit_status(145) .and. n.ge.47 .and. cyal.ge.2 .and. cyar.ge.2) then
        if (debugg) print *, 'screen: 145 oxazolidinone peptide'
        if (debugv) print *, 'vectar: 145 CCN1CC(OC1=O)C(O)C(CC2CCCCC2)NC=O'
        vectar(0) = 'CCN1CC(OC1=O)C(O)C(CC2CCCCC2)NC=O'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Oxazolidinone peptide',' ',oxazol_peptide)
            dummy = total(oxazol_peptide)
         endif
      endif
c
c SMARTS: styrene
c
      if (r187 .and. (audit_status(146) .or. audit_status(150)) 
     *  .and. n.ge.11 .and. cyar.ge.1 .and. at_st.ge.1. .and. db_cc.ge.1) then
        if (debugg) print *, 'screen: 150 styrene test'
        if (debugv) print *, 'vectar: 150 cC=CC'
        vectar(0) = 'cC=CC'
        call search( 0, 'ANYONE', sok )
        if (nhits .eq. 0) goto 1200
        if (debug) print *, 'Sty 0 '
c 
        if (r186 .and. audit_status(146) .and. 0.ne.index(ssmile,'#')) then
           if (debugg) print *, 'screen: 146 styrene test'
           if (debugv) print *, 'vectar: 146 C=CC#N'
           vectar(0) = 'C=CC#N'
           call search( 0, 'ANYONE', sok )
        else
           nhits = 0
        endif
c if styrene/cyano
        if (nhits .gt. 0) then
        if (debugg) print *, 'screen: 146 styrene/cyano'
        if (debug) print *, 'Sty 1 '
           if (debugv) print *, 'vectar: 146 c[C;H1;R0]=[C;H0;R0]C#N'
           vectar(0) = 'c[C;H1;R0]=[C;H0;R0]C#N'
           call search( 0, 'ANYONE', sok )
           if (debug) print *, 'Sty x', sok, nhits
           if (debug) print *, 'Sty y', nhits
           if (nhits .gt. 0) then
           if (debug) print *, 'Sty 1 - yes'
               call outlin(cclass,'SMARTS',' Styrene (di-substituted)',' ',styrene_mono)
               dummy = total(styrene_mono)
               if (debugv) print *, 'vectar: 146 [C;H1;R0]=[C;H0;R0](C#N)(C#N)'
               vectar(0) = '[C;H1;R0]=[C;H0;R0](C#N)(C#N)'
               call search( 0, 'ANYONE', sok )
               if (nhits .ne. 0) then
                  if (debug) print *, 'Sty 2 '
                  call outlin(cclass,'SMARTS',' Styrene (di-cyano)',' ',styrene_di)
                  dummy = total(styrene_di)
              endif
           endif
c styrene / cinnamamide
        elseif (r186 .and. audit_status(150) .and. al_carb.ge.1 .and. al_n.ge.1) then
        if (debug) print *, 'Sty 3 '
        if (debugg) print *, 'screen: 150 styrene/cinnamamide'
           if (debugv) print *, 'vectar: 150 c1ccccc1C=[C;R0]C(=O)N A'
           vectar(0) = 'c1ccccc1C=[C;R0]C(=O)N'
           call search( 0, 'ANYONE', sok )
           if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Styrene (cinnamamide)',' ',styrene_cinn)
              dummy = total(styrene_cinn)
              goto 1200
           else
              if (debugv) print *, 'vectar: 150 C(=O)C(=Cc1ccco1)N(=O)=O B'
              vectar(0) = 'C(=O)C(=Cc1ccco1)N(=O)=O'
              call search( 0, 'ANYONE', sok )
              if (nhits .ne. 0) then
                 call outlin(cclass,'SMARTS',' Styrene (di-substituted)',' ',0.2+styrene_mono)
                 dummy = total(0.2+styrene_mono)
              endif
           endif
        endif
 1200 continue
      endif
c
c SMARTS: piperazine/phenyl
c 
      if (r187 .and. audit_status(149) .and. piperazine .and. ar_cn.ge.1) then
        if (debugg) print *, 'screen: 149 piperazine/phenyl'
        if (debugv) print *, 'vectar: 149 [N;R1]1CC[N;R1](CC1)c2c[c;R1][c;R1][c;R1]c2'
        vectar(0) = '[N;R1]1CC[N;R1](CC1)c2c[c;R1][c;R1][c;R1]c2'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Piperazine/phenyl',' ', piper_phenyl)
            dummy = total(piper_phenyl)
         endif
      endif
c
c SMARTS: Glucopyranoside
c 
      if (r187 .and. audit_status(151) .and. al_o.ge.4 .and. at_oh.ge.2 .and. ar_eo.ge.1 .and. al_or.ge.1) then
        if (debugg) print *, 'screen: 151 glucopyranoside'
        if (debugv) print *, 'vectar: 151 [O;H]CC1OC(Oc2ccccc2)CCC1[O;H]'
        vectar(0) = '[O;H]CC1OC(Oc2ccccc2)CCC1[O;H]'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Glucopyranoside',' ', glucopyranoside)
            dummy = total(glucopyranoside)
         endif
      endif
c
c SMARTS: extended aromaticity
c 
      if (r187 .and. audit_status(152) .and. ar_carb.ge.2 .and. ar_n.ge.2 .and. ar_nh.ge.1 .and. 0.ne.index(ssmile,'c(=C')) then
      if (debugg) print *, 'screen: 152 extended aromaticity'
        if (debugv) print *, 'vectar: 152 C=c1nc(=O)c(=C)nc1=O'
        vectar(0) = 'C=c1nc(=O)c(=C)nc1=O'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Extended aromaticity',' ', extend_arom)
            dummy = total(extend_arom)
         endif
      endif
c
c SMARTS: pyrazine ortho
c 
      if (.not.r187 .or. .not. audit_status(164) .or. .not.pyrazine) goto 1300
      if (debugg) print *, 'screen: 164 pyrazine ortho'
      if (0.ne.index(ssmile,'#N') .or. 0.ne.index(ssmile,'N#')) then
         if (debugg) print *, 'screen: 164a pyrazine ortho / cyano'
         if (debugv) print *, 'vectar: 164 N#Cc1cnccn1'
         vectar(0) = 'N#Cc1cnccn1'
         call search( 0, 'UNIQUE', sok )
         if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Pyrazine ortho (cyano)',' ',nhits*pyr_ortho)
            dummy = total(nhits*pyr_ortho)
            goto 1300
         endif
      endif
      if (r187 .and. al_carb.ge.1 .and. al_o.ge.2) then
         if (debugg) print *, 'screen: 164b pyrazine ortho / carboxy'
         if (debugv) print *, 'vectar: 164 C(=O)(O)c1cnccn1'
         vectar(0) = 'C(=O)(O)c1cnccn1'
         call search( 0, 'UNIQUE', sok )
         if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Pyrazine ortho (carboxy)',' ',nhits*pyr_ortho*3/4)
            dummy = total(nhits*pyr_ortho*3/4)
            goto 1300
         endif
      endif
      if (r187 .and. at_halo.ge.1 .and. 0.eq.index(ssmile,'N(=O)=O') .and. 0.eq.index(ssmile,'O=N(=O)')) then
         if (debugg) print *, 'screen: 164c pyrazine ortho / halogen'
         if (debugv) print *, 'vectar: 164 [Cl,F,Br]c1cnccn1'
         vectar(0) = '[Cl,F,Br]c1cnccn1'
         call search( 0, 'UNIQUE', sok )
         if (nhits.ne.0) then
            call outlin(cclass,'SMARTS',' Pyrazine ortho (halogen)',' ',nhits*pyr_ortho/2)
            dummy = total(nhits*pyr_ortho/2)
            goto 1300
         endif
      endif
      if (r187 .and. al_s.ge.1) then
         if (debugg) print *, 'screen: 164d pyrazine ortho / thio'
         if (debugv) print *, 'vectar: 164 Sc1cnccn1'
         vectar(0) = 'Sc1cnccn1'
         call search( 0, 'UNIQUE', sok )
         if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Pyrazine ortho (thio)',' ',nhits*pyr_ortho/2)
            dummy = total(nhits*pyr_ortho/2)
            goto 1300
         endif
      endif
      if (r187 .and. ar_eo.ge.1 .and. al_c.ge.1) then
         if (debugg) print *, 'screen: 164e pyrazine ortho / ether'
         if (debugv) print *, 'vectar: 164 COc1cnccn1'
         vectar(0) = 'COc1cnccn1'
         call search( 0, 'UNIQUE', sok )
         if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Pyrazine ortho (ether)',' ',nhits*pyr_ortho/6)
            dummy = total(nhits*pyr_ortho/6)
         endif
      endif
 1300 continue
c
c SMARTS: gambogic acid analogs
c
      if (r187 .and. audit_status(165) .and. ncycle.ge.5 .and. al_carb.ge.1 .and. db_ccr.ge.1 .and. al_eo.ge.1 .and. fatom.ge.2 .and. cyal.ge.4 .and. cyar.ge.1) then
        if (debugg) print *, 'screen: 165 gambogic acid'
        if (debugv) print *, 'vectar: 165 C12CCOC1C(=O)CC=C2'
        vectar(0) = 'C12CCOC1C(=O)CC=C2'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
           call outlin(cclass,'SMARTS',' Gambogic acid analogs',' ',gambogic_acid)
           dummy = total(gambogic_acid)
        endif
      endif
c
c SMARTS: crown ethers
c
      if (r187 .and. audit_status(169) .and. cymax.ge.12) then
        if (debugg) print *, 'screen: 169 crown ethers'
        if (cyar .ge. 1) then
           if (debugv) print *, 'vectar: 169 OCCOc1ccccc1OCCOCCO A'
           vectar(0) = 'OCCOc1ccccc1OCCOCCO'
           call search( 0, 'ANYONE', sok )
           if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Crown ether (benzo)',' ',crown_benzo)
              dummy = total(crown_benzo)
           endif
        else 
           if (debugv) print *, 'vectar: 169 CCOCCOCCOCCO B'
           vectar(0) = 'CCOCCOCCOCCO'
           call search( 0, 'ANYONE', sok )
           if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Crown ether (simple)',' ',crown_ether)
              dummy = total(crown_ether)
           endif
        endif
      endif
c
c SMARTS: Di-aryl amide, di-ortho halogen
c
      if (r187 .and. audit_status(172) .and. ar_amide.ge.1 .and. (at_cl+at_f).ge.2 .and. cyar6.ge.2 .and. ar_n.ge.1) then
        if (debugg) print *, 'screen: 172 R7090'
        if (debugv) print *, 'vectar: 172 c1([Cl,F])cccc([Cl,F])c1C(=O)Nc'
        vectar(0) = 'c1([Cl,F])cccc([Cl,F])c1C(=O)Nc'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Di-aryl amide, di-ortho halogen',' ',0.80)
              dummy = total(0.80)
        endif
      endif
c
c SMARTS: Atropisomers (counter one 'hetero-aromatic extension' due to ortho twist [cf R7137])
c
      if (r187 .and. audit_status(173) .and. ar_n.ge.3 .and. cyar6.ge.2 .and. cyar5.ge.1 .and. arar .and. 0.ne.index(ssmile,'nn')) then
        if (debugg) print *, 'screen: 173 R7137'
        if (debugv) print *, 'vectar: 173 c1ccc(cc1)[n;R1]2c[n;H0][n;H0]c2c3ccc[n,c]c3'
        vectar(0) = 'c1ccc(cc1)[n;R1]2c[n;H0][n;H0]c2c3ccc[n,c]c3'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Atropisomers',' ',-0.31)
              dummy = total(-0.31)
        endif
      endif
c
c SMARTS: Methaqualone analogs
c
      if (r187 .and. audit_status(174) .and. cyar.ge.3 .and. cy6.ge.3 .and. ar_carb.ge.1 .and. ar_n.ge.1 .and. fatom.ge.2 .and. (at_cl.ge.1 .or. at_bz.ge.1)) then
        if (debugg) print *, 'screen: 174 methaqualone analogs'
        if (debugv) print *, 'vectar: 174 [C,Cl;R0]c1ccccc1[n;R1]3c[c,n]ccc3=O'
        vectar(0) = '[C,Cl;R0]c1ccccc1[n;R1]3c[c,n]ccc3=O'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Methaqualone analogs',' ',-0.60)
              dummy = total(-0.60)
        endif
      endif
c
c SMARTS: Nitrogen bridge across two 5-rings
c
      if (r187 .and. audit_status(176) .and. al_n.ge.2 .and. cyal5.ge.2 .and. fatom.ge.2) then
        if (debugg) print *, 'screen: 176 nitrogen bridge across two 5-rings'
        if (debugv) print *, 'vectar: 176 C1NCC2CNCC12'
        vectar(0) = 'C1NCC2CNCC12'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' nitrogen bridge over 5-rings',' ',0.75)
              dummy = total(0.75)
        endif
      endif
c
c SMARTS: Iodide ortho to in-ring nitrogen (on 6-ring)
c
      if (r187 .and. audit_status(178) .and. ar_n.ge.1 .and. cy6.ge.1 .and. 0.ne.index(ssmile,'I')) then
        if (debugg) print *, 'screen: 178 iodide ortho to in-ring nitrogen (on 6-ring)'
        if (debugv) print *, 'vectar: 178 Ic1[n;R1][c;R1]ccc1'
        vectar(0) = 'Ic1[n;R1][c;R1]ccc1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' I/n ortho offset',' ',0.30)
              dummy = total(0.30)
        endif
      endif
c
c SMARTS: R7447 scaffold needs YCCCY, but use this instead
c
      if (r187 .and. audit_status(180) .and. al_eo.ge.1 .and. cyal6.ge.1 .and. al_amide.ge.1. and. al_n.ge.3 
     *   .and. 0.ne.index(ssmile,'C#N') ) then
        if (debugg) print *, 'screen: 180 R7447 substitute for YCCCY'
        if (debugv) print *, 'vectar: 180 N#CCNC(=O)C1(N)CCOCC1'
        vectar(0) = 'N#CCNC(=O)C1(N)CCOCC1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' para-di-fragment pyran',' ',1.00)
              dummy = total(1.00)
        endif
      endif

c
c SMARTS: R7583 with pyrrolidine-phenyl ortho twist
c
      if (r187 .and. audit_status(182) .and. pyrazole .and. arar .and. cyar6.ge.1 .and. cyar5.ge.1) then
        if (at_bz.ge.1 .and. ch3.ge.1) then
           if (debugg) print *, 'screen: 182 pyrrolidine-phenyl ortho twist A'
           if (debugv) print *, 'vectar: 182 c1cccc[c;R1]1n2c([C;H3])ccn2 A'
           vectar(0) = 'c1cccc[c;R1]1n2c([C;H3])ccn2'
           call search( 0, 'ANYONE', sok )
           if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Pyrrolidine-phenyl ortho #1',' ',-.80)
              dummy = total(-.80)
              goto 1500
           endif
        endif
        if (at_cl.ge.1) then
           if (debugg) print *, 'screen: 182 pyrrolidine-phenyl ortho twist B'
           if (debugv) print *, 'vectar: 182 c1cccc(Cl)[c;R1]1n2nccc2 B'
           vectar(0) = 'c1cccc(Cl)[c;R1]1n2nccc2'
           call search( 0, 'ANYONE', sok )
           if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Pyrrolidine-phenyl ortho #2',' ',-.50)
              dummy = total(-.50)
           endif
        endif
      endif
 1500 continue

c
c SMARTS: R7592 multiple meth-oxy
c   Better:    24     gain     6.53   avg gain  0.27
c   Worse:    125     loss    64.16   avg loss  0.51
c   Totals:   149      net   -57.63
c
c     if (r187 .and. audit_status(xxx) .and. ch3.ge.2 .and. cy6.ge.1 .and. cyar.ge.1 .and. 
c    *     (0.ne.index(ssmile,'COc') .or. 0.ne.index(ssmile,'c(OC)')) ) then
c       vectar(0) = 'c1(O[C;H3])c(O[C;H3])cccc1'
c       call search( 0, 'ANYONE', sok )
c       if (nhits .ne. 0) then
c          call outlin(cclass,'SMARTS',' Meth-oxy pair',' ',-.55)
c          dummy = total(-.55)
c          vectar(0) = 'c1(O[C;H3])c(O[C;H3])c(O[C;H3])ccc1'
c          call search( 0, 'ANYONE', sok )
c          if (nhits .ne. 0) then
c             call outlin(cclass,'SMARTS',' Meth-oxy pair (additional)',' ',-.15)
c             dummy = total(-.15)
c          endif
c       endif
c     endif

c
c SMARTS: ortho-phenols
c
      if (r187 .and. audit_status(99) .and. arar .and. at_oh.ge.1 .and. ar_o.ge.1 .and. cyar5.ge.1 .and. cyar6.ge.1) then 
        if (debugg) print *, 'screen: 99 ortho-phenols'
        if (debugv) print *, 'vectar: 99 c1cccc([O;H])c1-co'
        vectar(0) = 'c1cccc([O;H])c1-co'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' ortho-phenols',' ',0.90)
              dummy = total(0.90)
        endif
      endif

c
c SMARTS: CF3 on diazo-oxy
c
      if (r187 .and. audit_status(107) .and. at_f.ge.3 .and. ar_o.ge.1 .and. cy5.ge.1) then 
        if (debugg) print *, 'screen: 107 CF3 on diazo-oxy'
        if (debugv) print *, 'vectar: 107 C(F)(F)(F)co'
        vectar(0) = 'C(F)(F)(F)co'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' CF3 on diazo-oxy',' ',0.85)
              dummy = total(0.85)
        endif
      endif

c
c SMARTS: Podophyllotoxins : TOPOsides
c
      if (r187 .and. audit_status(108) .and. ncycle.ge.4 .and. cyal.ge.3 .and. al_or.ge.3 .and. cy5.ge.2 .and. cy6.ge.2 .and. cyar.ge.1) then 
        if (debugg) print *, 'screen: 108 podophyllotoxins'
        if (debugv) print *, 'vectar: 108 C4Oc3cc2CC1COCC1Cc2cc3O4'
        vectar(0) = 'C4Oc3cc2CC1COCC1Cc2cc3O4'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Podophyllotoxins/TOPOsides',' ',1.70)
              dummy = total(1.70)
        endif
      endif

c
c SMARTS: di-methoxy offsetting correction
c
      if (r187 .and. audit_status(109) .and. ar_eo.ge.2 .and. ch3.ge.2) then 
        if (debugg) print *, 'screen: 109 di-methoxy'
        if (debugv) print *, 'vectar: 109a [c;H]c(O[CH3])c(O[CH3])[c;H]'
        vectar(0) = '[c;H1]c(O[C;H3])c(O[C;H3])[c;H1]'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' di-methoxy',' ',-0.20)
              dummy = total(-0.20)
        endif
      endif

c
c SMARTS: tri-methoxy offsetting correction
c
      if (r187 .and. audit_status(109) .and. ar_eo.ge.3 .and. ch3.ge.3) then 
        if (debugg) print *, 'screen: 109 tri-methoxy'
        if (debugv) print *, 'vectar: 109b c(O[CH3])c(O[CH3])c(O[CH3])'
        vectar(0) = 'c1c(O[C;H3])c(O[C;H3])c(O[C;H3])cc1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' tri-methoxy',' ',0.10)
              dummy = total(0.10)
        endif
      endif

c
c SMARTS: Quat-carbon with hydroxy, in di-hetero ring environment
c
      if (r187 .and. audit_status(136) .and. ar_s.ge.2 .and. cyar.ge.2) then 
        if (debugg) print *, 'screen: 136 quat-carbon'
        if (debugv) print *, 'vectar: 136 CC(O)(c1nccs1)c2nncs2'
        vectar(0) = 'CC(O)(c1nccs1)c2nncs2'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' quat-carbon, OH/di-hetero',' ',2.00)
              dummy = total(2.00)
        endif
      endif
c
c SMARTS: para-hydroxy piperidine
c
      if (r187 .and. audit_status(153) .and. piperidine .and. at_oh.ge.1 .and. nfrag.ge.2) then
        if (debugg) print *, 'screen: 153 para-hydroxy piperidine'
        if (debugv) print *, 'vectar: 153 [O;H][C;H1]1[C;H2][C;H2]N(C=O)[C;H2][C;H2]1'
        vectar(0) = '[O;H][C;H1]1[C;H2][C;H2]N(C=O)[C;H2][C;H2]1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' para-hydroxy piperidine',' ',0.80)
              dummy = total(0.80)
        endif
      endif
c
c SMARTS: conjugated triple bonds
c
      if (r187 .and. audit_status(155) .and. 0.ne.index(ssmile,'C#CC#CC#C')) then
        if (debugg) print *, 'screen: 155 conjugated triple bonds'
        if (debugv) print *, 'vectar: 155 C#CC#CC#C'
        vectar(0) = 'C#CC#CC#C'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' conjugated triple bonds',' ',1.00)
              dummy = total(1.00)
        endif
      endif
c
c SMARTS: Scaffold: offset CF3 Y-X proximity in multi-Y case
c
      if (r187 .and. audit_status(166) .and. at_cf3.ge.1 .and. nfrag.ge.5) then
        if (debugg) print *, 'screen: 166 CF3 multi-Y offset'
        if (debugv) print *, 'vectar: 166 FC(F)(F)C([N,O,S])[C;H2][N,O,S]'
        vectar(0) = 'FC(F)(F)C([N,O,S])[C;H2][N,O,S]'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' CF3 multi-Y offset',' ',-1.0)
              dummy = total(-1.0)
        endif
      endif
c
c SMARTS: Scaffold: Sulfonyl benzyl to hetero ring
c
      if (r187 .and. audit_status(184) .and. nfrag.ge.3 .and. at_bz.ge.1 .and. ch3.ge.1 .and. s184.ge.1) then
        if (debugg) print *, 'screen: 184 Sulfonyl benzyl to hetero ring'
        if (debugv) print *, 'vectar: 184 S(=O)(=O)Ccn'
        vectar(0) = 'S(=O)(=O)Ccn'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Sulfonyl benzyl to hetero ring',' ',0.7)
              dummy = total(0.7)
        endif
      endif
c
c SMARTS: Scaffold: Amide on indole
c
      if (r187 .and. audit_status(189) .and. ar_amide.ge.1 .and. ar_nh.ge.1 .and. fatom.ge.2 .and. cyar6.ge.1 .and. cyar6.ge.1) then
        if (debugg) print *, 'screen: 189 Amide on indole'
        if (debugv) print *, 'vectar: 189 NC(=O)c2cc1ccccc1[n;H]2'
        vectar(0) = 'NC(=O)c2cc1ccccc1[n;H]2'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' Amide on indole',' ',0.37)
              dummy = total(0.37)
        endif
      endif
c
c SMARTS: Scaffold: R7949
c
      if (r187 .and. audit_status(190) .and. cy5.ge.2 .and. ar_o.ge.1 .and. ar_n.ge.1 .and. cyal.ge.1) then
        if (debugg) print *, 'screen: 190 R7949'
        if (debugv) print *, 'vectar: 190 c1[n,o]nc(NC(=O)[N,O])c1'
        vectar(0) = 'c1[n,o]nc(NC(=O)[N,O])c1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' R7949',' ',1.20)
              dummy = total(1.20)
        endif
      endif
c
c SMARTS: Scaffold: R7962 hydrate
c
      if (r187 .and. audit_status(191) .and. 0.ne.index(ssmile,'C#N') .and. cy6.ge.2 .and. ar_n.ge.2) then
        if (debugg) print *, 'screen: 191 R7962'
        if (debugv) print *, 'vectar: 191 c1nc(cnc1)c2ccccc2C#N'
        vectar(0) = 'c1nc(cnc1)c2ccccc2C#N'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' R7962 hydrate',' ',0.70)
              dummy = total(0.70)
        endif
      endif
c
c SMARTS: Scaffold: R7581 H-bond potential
c
      if (r187 .and. audit_status(194) .and. cyar5.ge.1 .and. ar_s.ge.1 .and. ar_n.ge.1) then
        if (debugg) print *, 'screen: 194 R7581'
        if (debugv) print *, 'vectar: 194 NC(=O)c1ncsc1C(O)=O'
        vectar(0) = 'NC(=O)c1ncsc1C(O)=O'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' R7581 H-bond potential',' ',2.00)
              dummy = total(2.00)
        endif
      endif
c
c SMARTS: Scaffold: R8118 H-bond potential
c
      if (r187 .and. audit_status(195) .and. cyar5.ge.1 .and. ar_n.ge.1 .and. ar_o.ge.1 .and. al_carb.ge.1 .and. al_o.ge.1 .and. al_n.ge.1) then
        if (debugg) print *, 'screen: 195 R8118'
        if (debugv) print *, 'vectar: 195 onc([O;H])c[C;H2][C;H1][N;H2]'
        vectar(0) = 'onc([O;H])c[C;H2][C;H1][N;H2]'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' R8118 H-bond potential',' ',1.00)
              dummy = total(1.00)
        endif
      endif
c
c SMARTS: Scaffold: 1,3,4 oxadiazole - phenyl
c
      if (r187 .and. audit_status(198) .and. cyar5.ge.1 .and. cyar6.ge.1 .and. ar_n.ge.2 .and. ar_o.ge.1) then
        if (debugg) print *, 'screen: 198 1,3,4 oxadiazole - phenyl'
        if (debugv) print *, 'vectar: 198 n1ncoc1-[c;R1]'
        vectar(0) = 'n1ncoc1-[c;R1]'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' 1,3,4 oxadiazole - phenyl',' ',0.60)
              dummy = total(0.60)
        endif
      endif
c
c SMARTS: Scaffold: oxazole
c
      if (r187 .and. audit_status(199) .and. cyar5.ge.1 .and. ar_n.ge.1 .and. ar_o.ge.1) then
        if (debugg) print *, 'screen: 199 oxazole'
        if (debugv) print *, 'vectar: 199 n1co[c;R][c;R1]1'
        vectar(0) = 'n1co[c;R][c;R1]1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' oxazole',' ',0.30)
              dummy = total(0.30)
        endif
      endif
c
c SMARTS: Scaffold: ethers ortho to amide
c
      if (r187 .and. audit_status(200) .and. ar_eo.ge.2 .and. cyar6.ge.1 .and. al_o.ge.3 .and. al_n.ge.1) then
        if (debugg) print *, 'screen: 200 ethers ortho to amide'
        if (debugv) print *, 'vectar: 200 O(C)c1cccc(OC])c1C([N;H])=O'
        vectar(0) = 'O(C)c1cccc(OC)c1C([N;H])=O'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
              call outlin(cclass,'SMARTS',' ethers ortho to amide',' ',-0.40)
              dummy = total(-0.40)
        endif
      endif
c
c SMARTS: Scaffold: methyoxy ortho pair on benzoquinone 
c
      if (r187 .and. audit_status(201) .and. al_eo.ge.2 .and. al_carb.ge.2 .and. cyal6.ge.1 .and. db_ccr.ge.2) then
        if (debugg) print *, 'screen: 201 benzoquinone methyoxy ortho pair'
        if (debugv) print *, 'vectar: 201 O=C1[C;R1]=CC(=O)C(O[C;H3])=C1O[C;H3]'
        vectar(0) = 'O=C1[C;R1]=CC(=O)C(O[C;H3])=C1O[C;H3]'
        call search( 0, 'UNIQUE', sok )
        if (nhits .eq. 2) then
              call outlin(cclass,'SMARTS',' benzoquinone methyoxy ortho pair',' ',0.70)
              dummy = total(0.70)
        elseif (nhits .eq. 4) then
              call outlin(cclass,'SMARTS',' benzoquinone methyoxy ortho pair',' ',2*0.70)
              dummy = total(2*0.70)
        endif
      endif

c=================================================================================
c END OF SMARTS
c=================================================================================

c
c some de-activated corrections
c

c
c SMARTS: vinyl conjugation
c
c     if (r187 .and. audit_status(xx) .and. 0.ne.index(ssmile,'C=C')) then
c       vectar(0) = 'O=C[C;H1]=[C;R0;H1]c1ccccc1'
c       call search( 0, 'ANYONE', sok )
c       if (nhits .ne. 0) then
c             call outlin(cclass,'SMARTS',' vinyl conjugation',' ',0.20)
c             dummy = total(0.20)
c       endif
c     endif

c
c SMARTS: allow H-bond between SOH and amine
c
c     if (audit_status(xx) .and. cyar.ge.1 .and. 0.ne.index(ssmile,'N') .and. 0.ne.index(ssmile,'S')) then
c       vectar(0) = 'S([O;H])cc[N;H2]'
c       call search( 0, 'UNIQUE', sok )
c       if (nhits .ne. 0) then
c       if (1.eq.1 .or. debugg) print *, 'screen: xx', nhits
c             call outlin(cclass,'SMARTS',' sulfonic H-bond',' ',nhits*0.50)
c             dummy = total(nhits*0.50)
c       endif
c     endif

c
c SMARTS: extra benzyl
c
c     if (audit_status(xxx)) then
c       vectar(0) = 'c1cccc([O,N])c1[C;R0](C)(C)'
c       call search( 0, 'ANYONE', sok )
c       if (nhits .ne. 0) then
c           call outlin(cclass,'SMARTS',' extra benzyl',' ',-0.10)
c           dummy = total(-0.10)
c        endif
c     endif
c

c=================================================================================
c Tests that follow are just warnings
c=================================================================================
c
c SMARTS: erythromycin
c
      if (r187 .and. audit_status(52) .and. cyal.ge.3 .and. cymax.ge.14) then
        if (debugg) print *, 'screen: 52 erythromycin warning'
        if (debugv) print *, 'vectar:  52 C2CCCOC2OC1CC(=O)OCCCCC(=O)CCCC(OC3OCCCC3)C1'
        vectar(0) = 'C2CCCOC2OC1CC(=O)OCCCCC(=O)CCCC(OC3OCCCC3)C1'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Erythromycin -- may be low','Warning',0.015)
            dummy = total(0.015)
         endif
      endif
c
c SMARTS: polypeptide warning
c
      if (r187 .and. audit_status(120) .and. n.ge.18 .and. al_amide.ge.3) then
        if (debugg) print *, 'screen: 120 pp flexible warning'
        vectar(0) = 'NC[C;R0](=O)NC[C;R0](=O)NC[C;R0](=O)'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Polypeptides unpredictable (flexible)','Warning',0.015)
            dummy = total(0.015)
         endif
      endif
c
c SMARTS: polypeptide/macrocycle warning
c
      if (r187 .and. audit_status(121) .and. cyal.ge.1. .and. cymax.gt.14) then
        if (debugg) print *, 'screen: 121 macrocycle warning'
        vectar(0) = 'NC[C;R1](=O)NC[C;R1](=O)NC[C;R1](=O)NC[C;R1](=O)'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Peptide macrocycles unpredictable','Warning',0.015)
            dummy = total(0.015)
        endif
      endif
c
c SMARTS: allopurinol/guanine tautomers warning
c
      if (r187 .and. audit_status(171) .and. cyar.ge.2 .and. ar_carb.ge.1 .and. ar_n.ge.4 .and. fatom.ge.2 .and. cy5.ge.1 .and. cy6.ge.1) then
        if (debugg) print *, 'screen: 171 tautomer warning'
        vectar(0) = 'O=c1[n]cnc2[n](C)ncc12'
        call search( 0, 'ANYONE', sok )
        if (nhits .ne. 0) then
            call outlin(cclass,'SMARTS',' Tautomeric form may be higher','Keto/enol',0.000)
            dummy = total(0.000)
        else
            vectar(0) = 'Nc2nc1[n]cnc1c(=O)[n]2 '
            call search( 0, 'ANYONE', sok )
            if (nhits .ne. 0) then
                call outlin(cclass,'SMARTS',' Tautomeric form may be higher','Keto/enol',0.000)
                dummy = total(0.000)
            endif
        endif
      endif

c
c warn of very low/high estimated values (error levels 51 and 52)
c
      dummy = total( 0.0 )
      clogp_max = 12.0
      clogp_min = -1 * clogp_max
      if (audit_status(170) .and. dummy.gt.clogp_max) then
         ccap = clogp_max - dummy
         dsink = total(ccap)
         call outlin(cclass,ctype,' ClogP exceeds maximum, capped ','Warning',ccap)
         call outlin(cclass,ctype,' Not predictive for transport or binding','Warning',0.0)
      endif

      if (audit_status(170) .and. dummy.lt.clogp_min) then
         ccap = -1 * (dummy - clogp_min)
         dsink = total(ccap)
         call outlin(cclass,ctype,' ClogP falls below minimum, capped ','Warning',ccap)
         call outlin(cclass,ctype,' Not predictive for transport or binding','Warning',0.0)
      endif

      if (dummy .gt. 8.0 ) then
         iflag = 51
         if (warn) then
            call outlin(cclass,ctype,' Value > 8 not applicable in modeling','Warning',0.0)
         endif
      endif
      if (dummy .lt. -3.5 ) then
         iflag = 52
         if (warn) then
            call outlin(cclass,ctype,' Value < -3.5 not applicable in modeling','Warning',0.0)
         endif
      endif
c
      return
      end
