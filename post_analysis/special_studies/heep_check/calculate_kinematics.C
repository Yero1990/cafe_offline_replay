void calculate_kinematics(){

  /* Brief: this code defines the calculated (e,e'p) elastics kinematics
     using minimum required initial information, e.g.
     (Ebeam, e- momentum) or (Ebeam, e- angle) or (e- momentum, e- angle)

     The quantities are:
     kf_calc (Eb, the)
     Pf_calc (Eb. the)
     thp_calc (Eb, the)
     
     Then the calculated-measured quantities are filled into histograms for fitting, e.g.,
     dkf = (kf_calc - kf_meas),
     dPf =  (Pf_calc - Pf_meas),
     dthp = (thp_calc - thp_meas)

     keep in mind that the measured quantities do have energy loss whereas
     the calculated quantities do not, and eloss will have to be subtarcted
     from the mesured quantities. To achieve this, calculated the difference between
     the reconstructed and thrown variables from SIMC, then a fit to that difference
     should give the numerical value of eloss

   */

}
