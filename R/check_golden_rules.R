check_golden_rules<-function(curformula,NOPS_check=FALSE){


  numnitrogens<-check_element(curformula,"N")
  numcarbons<-check_element(curformula,"C")

  numoxygens<-check_element(curformula,"O")

  numhydrogens<-check_element(curformula,"H")

  numphos<-check_element(curformula,"P")
  numsulphur<-check_element(curformula,"S")



  if(isTRUE((numcarbons<1)[1])){

    bool_check<-0
  }else{

    max_hydrogens<-(2*numcarbons)+numnitrogens+2

    nitrogen_to_carbon_ratio<-numnitrogens/numcarbons

    oxygen_to_carbon_ratio<-numoxygens/numcarbons

    phosphorus_to_carbon_ratio<-numphos/numcarbons

    sulphur_to_carbon_ratio<-numsulphur/numcarbons

    hydrogens_to_carbon_ratio<-numhydrogens/numcarbons

    bool_check<-1



    if(isTRUE((hydrogens_to_carbon_ratio<0.1 | hydrogens_to_carbon_ratio>6)[1])){

      bool_check=0
    }

    if(isTRUE((nitrogen_to_carbon_ratio>4)[1])){

      bool_check=0
    }

    if(isTRUE((oxygen_to_carbon_ratio>3)[1])){

      bool_check=0
    }

    if(isTRUE((phosphorus_to_carbon_ratio>2)[1])){

      bool_check=0
    }

    if(isTRUE((sulphur_to_carbon_ratio>3)[1])){

      bool_check=0
    }

    if(isTRUE((NOPS_check==TRUE)[1])){
      #NOPS>1
      if(isTRUE((numnitrogens>1 & numoxygens>1 & numphos>1 & numsulphur>1)[1])){
        if(isTRUE((numnitrogens>10 | numoxygens>20 | numphos>4 | numsulphur>3)[1])){

          bool_check<-0
        }

      }


      #NOP>3
      if(isTRUE((numnitrogens>3 & numoxygens>3 & numphos>3)[1])){
        if(isTRUE((numnitrogens>11 | numoxygens>22 | numphos>6)[1])){

          bool_check<-0
        }

      }

      #OPS>1
      if(isTRUE((numoxygens>1 & numphos>1 & numsulphur>1)[1])){
        if(isTRUE((numoxygens>14 | numphos>3 | numsulphur>3)[1])){

          bool_check<-0
        }

      }

      #PSN>1
      if(isTRUE((numnitrogens>1 & numphos>1 & numsulphur>1)[1])){
        if(isTRUE((numnitrogens>4 | numphos>3 | numsulphur>3)[1])){

          bool_check<-0
        }

      }

      #NOS>6
      if(isTRUE((numnitrogens>6 & numoxygens>6 & numsulphur>6)[1])){
        if(isTRUE((numnitrogens>19 | numoxygens>14 | numsulphur>8)[1])){

          bool_check<-0
        }

      }

    }


  }
  res<-cbind(curformula,bool_check)
  res<-as.data.frame(res)


  return(res)

}
