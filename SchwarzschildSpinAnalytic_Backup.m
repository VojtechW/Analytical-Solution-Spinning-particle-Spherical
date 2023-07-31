(* ::Package:: *)

(* ::Chapter:: *)
(*Ancillary file to Analytical solutions for the motion of spinning particles near spherically symmetric black holes and exotic compact objects*)


(* ::Text:: *)
(*Created by Vojtěch Witzany and Gabriel Andres Piovano.*)


(* ::Text:: *)
(*This file contains formulas for the *)
(*- Mino frequencies (Υ^(t,r,ϑ,ψ)) ,*)
(*- spin-frequency corrections s∥δΥ at fixed coordinate semi-latus rectum p and eccentricity e,*)
(*- the explicit solutions for the motion x=t,r,ϑ,φ,ψ are given as x=x0 + s∥ δx in terms of the angle variables q^(t,r,φ,ψ) (defined with the help of an auxilliary function χ0(q^r)),*)
(*- explicit solution for the spin vector s*)
\[CapitalUpsilon]t0 = (Sqrt[((-4*e^2 + (-2 + p)^2)*p)/(-3 - e^2 + p)]*
      (-((p*(-6 + 2*e + p)*EllipticE[(4*e)/(-6 + 2*e + p)])/
         ((-1 + e^2)*(-4 + p))) + (p*(36 + 2*e^2*(-2 + p) + (-14 + p)*p)*
         EllipticK[(4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^2) + 
       2*(6 + 2*e - p)*(((-8 + e^2*(8 - 3*p) + (-1 + p)*p)*
           EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), 
            (4*e)/(-6 + 2*e + p)])/((-1 + e)*(1 + e)^2*(-4 + p)^2) + 
         EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
           (4*e)/(-6 + 2*e + p)]/(-2 - 2*e + p))))/
     (2*EllipticK[(4*e)/(-6 + 2*e + p)])
 
\[CapitalUpsilon]r0 = (Sqrt[-((p*(-6 + 2*e + p))/(3 + e^2 - p))]*Pi)/
     (2*EllipticK[(4*e)/(-6 + 2*e + p)])
 
\[CapitalUpsilon]\[CurlyPhi]0 = p/Sqrt[-3 - e^2 + p]
 
\[CapitalUpsilon]\[Psi]0 = 
    -(((-1 + e)*(6 + 2*e - p)*Im[EllipticPi[
          (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
            (1 + e + I*Sqrt[-3 - e^2 + p])), (4*e)/(-6 + 2*e + p)]] + 
       Sqrt[-3 - e^2 + p]*(4*EllipticK[(4*e)/(-6 + 2*e + p)] + 
         (-6 - 2*e + p)*Re[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + 
                  p]))/((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), 
            (4*e)/(-6 + 2*e + p)]]))/
      (Sqrt[((-4*e^2 + (-2 + p)^2)*(-3 - e^2 + p))/p]*
       EllipticK[(4*e)/(-6 + 2*e + p)]))
 
\[Delta]\[CapitalUpsilon]t = 
    ((4*(1 + e)*(3 + e^2 - p)*(-4 + p)^2*p*(-2 - 2*e + p)^2*(-6 + 2*e + p)*
        (-2 + 2*e + p)*EllipticE[(4*e)/(-6 + 2*e + p)]^2)/(6 + 2*e - p) - 
      4*(-4 + p)*(-2 - 2*e + p)*EllipticE[(4*e)/(-6 + 2*e + p)]*
       ((1 + e)*p*(-6 + 2*e + p)*(-28 + e^4*(4 + p) - 2*e^2*(-3 + p)*
           (4 + p) + p*(33 + 2*(-7 + p)*p))*EllipticK[(4*e)/(-6 + 2*e + p)] + 
        2*(3 + e^2 - p)*(-2 + 2*e + p)*((2 + 2*e - p)*(8 + p - p^2 + 
            e^2*(-8 + 3*p))*EllipticPi[(2*e*(-4 + p))/((1 + e)*
              (-6 + 2*e + p)), (4*e)/(-6 + 2*e + p)] + (-1 + e)*(1 + e)^2*
           (-4 + p)^2*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
            (4*e)/(-6 + 2*e + p)])) + EllipticK[(4*e)/(-6 + 2*e + p)]*
       (-((1 + e)*(2 + 2*e - p)*p*(-3072 + 48*e^6*p + 
           e^4*p*(112 + p*(32 + 3*(-12 + p)*p)) + 
           p*(4560 + p*(-2848 + p*(892 + 7*(-19 + p)*p))) + 
           e^2*(3072 + p*(-2672 + p*(1024 + p*(-280 + (50 - 3*p)*p)))))*
          EllipticK[(4*e)/(-6 + 2*e + p)]) + (6 + 2*e - p)*
         (8*(2 + 2*e - p)*(-128 + 3*e^6*p^2 + e^4*(-128 + p*(80 + 3*p - 
                6*p^2)) - p*(-144 + p*(55 + p*(2 + (-6 + p)*p))) + 
            e^2*(256 + p*(-224 + p*(113 - 40*p + 6*p^2))))*
           EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), 
            (4*e)/(-6 + 2*e + p)] + (-1 + e)*(1 + e)^2*(-4 + p)^3*
           (32 - p*(32 + (-11 + p)*p) + e^2*(-32 + p^2))*
           EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
            (4*e)/(-6 + 2*e + p)])))/(16*(-1 + e)*(1 + e)^2*(4 - p)^3*
      (-2 - 2*e + p)*(-3 - e^2 + p)^(3/2)*EllipticK[(4*e)/(-6 + 2*e + p)]^2)
 
\[Delta]\[CapitalUpsilon]r = 
    (Sqrt[(-4*e^2 + (-2 + p)^2)*(-6 + 2*e + p)*(-3 - e^2 + p)]*Pi*
      ((3 + e^2 - p)*EllipticE[(4*e)/(-6 + 2*e + p)] + 
       (-6 - 2*e + p)*EllipticK[(4*e)/(-6 + 2*e + p)]))/
     (4*(6 + 2*e - p)*(3 + e^2 - p)^2*EllipticK[(4*e)/(-6 + 2*e + p)]^2)
 
\[Delta]\[CapitalUpsilon]\[CurlyPhi] = 
    -1/2*((3 + e^2)*Sqrt[-4 + (4 - 4*e^2)/p + p])/(-3 - e^2 + p)^(3/2)
 
t0 = qt - (Sqrt[(-4*e^2 + (-2 + p)^2)/(-6 + 2*e + p)]*qr*
       (-((p*(-6 + 2*e + p)*EllipticE[(4*e)/(-6 + 2*e + p)])/
          ((-1 + e)*(1 + e)*(-4 + p))) + 
        (p*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*EllipticK[
           (4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^2) - 
        (2*(6 + 2*e - p)*(8 + p - p^2 + e^2*(-8 + 3*p))*
          EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), 
           (4*e)/(-6 + 2*e + p)])/((-1 + e)*(1 + e)^2*(-4 + p)^2) - 
        (2*(6 + 2*e - p)*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
           (4*e)/(-6 + 2*e + p)])/(2 + 2*e - p)))/Pi + 
     Sqrt[(-4*e^2 + (-2 + p)^2)/(-6 + 2*e + p)]*
      (-((p*(-6 + 2*e + p)*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)])/
         ((-1 + e)*(1 + e)*(-4 + p))) + (p*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*
         EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^2) - 
       (2*(6 + 2*e - p)*(8 + p - p^2 + e^2*(-8 + 3*p))*
         EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), \[Chi]0, 
          (4*e)/(-6 + 2*e + p)])/((-1 + e)*(1 + e)^2*(-4 + p)^2) - 
       (2*(6 + 2*e - p)*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
          \[Chi]0, (4*e)/(-6 + 2*e + p)])/(2 + 2*e - p) + 
       (e*p*Sqrt[(-6 + 2*e + p)*(-6 + p + 2*e*Cos[2*\[Chi]0])]*
         Sin[2*\[Chi]0])/((-1 + e^2)*(-6 + 2*e^2 + p + 
          e*(-4 + p)*Cos[2*\[Chi]0])))
 
r0 = (p*(-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
          (4*e)/(-6 + 2*e + p)]^2))/((1 + e)*(-6 + 2*e + p) - 
      2*e*(-4 + p)*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
         (4*e)/(-6 + 2*e + p)]^2)
 
\[CurlyPhi]0 = (Sqrt[p*(-6 + 2*e + p)]*Pi*q\[CurlyPhi] + 
      2*p*Pi*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)] - 
      2*p*qr*EllipticK[(4*e)/(-6 + 2*e + p)])/(Sqrt[p*(-6 + 2*e + p)]*Pi)
 
\[Psi]0 = (Sqrt[(-3 - e^2 + p)/(-6 + 2*e + p)]*
      ((2 + 2*e - p)*(-2 + 2*e + p)^2*Sqrt[(-6 + 2*e + p)/(-3 - e^2 + p)]*Pi*
        q\[Psi] + 2*Sqrt[-4*e^2 + (-2 + p)^2]*qr*
        (-4*(-2 + 2*e + p)*EllipticK[(4*e)/(-6 + 2*e + p)] - 
         ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*
           Im[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
              ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), 
             (4*e)/(-6 + 2*e + p)]])/Sqrt[-3 - e^2 + p] + 
         (4*(1 + e)^2 - (-4 + p)^2)*Re[EllipticPi[
            (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
              (1 + e + I*Sqrt[-3 - e^2 + p])), (4*e)/(-6 + 2*e + p)]]) - 
       2*Sqrt[-4*e^2 + (-2 + p)^2]*Pi*(-4*(-2 + 2*e + p)*EllipticF[\[Chi]0, 
           (4*e)/(-6 + 2*e + p)] - ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*
           Im[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
              ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), \[Chi]0, 
             (4*e)/(-6 + 2*e + p)]])/Sqrt[-3 - e^2 + p] + 
         (4*(1 + e)^2 - (-4 + p)^2)*Re[EllipticPi[
            (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
              (1 + e + I*Sqrt[-3 - e^2 + p])), \[Chi]0, 
            (4*e)/(-6 + 2*e + p)]])))/((4*e^2 - (-2 + p)^2)*(-2 + 2*e + p)*Pi)
 
\[Chi]0qr = JacobiAmplitude[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
     (4*e)/(-6 + 2*e + p)]
 
Sv0 = {(p*Sqrt[s^2 - sp^2]*(-6 + 2*e + p - 
        4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)]^2)*
       ((Cos[(Sqrt[(-3 - e^2 + p)/(-6 + 2*e + p)]*((2 + 2*e - p)*(-2 + 2*e + 
                 p)^2*Sqrt[(-6 + 2*e + p)/(-3 - e^2 + p)]*Pi*q\[Psi] + 
              2*Sqrt[-4*e^2 + (-2 + p)^2]*qr*(-4*(-2 + 2*e + p)*EllipticK[
                  (4*e)/(-6 + 2*e + p)] - ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + 
                   p)*Im[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
                     ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), 
                    (4*e)/(-6 + 2*e + p)]])/Sqrt[-3 - e^2 + p] + 
                (4*(1 + e)^2 - (-4 + p)^2)*Re[EllipticPi[(2*e*(-4 + p + 
                      (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*(1 + e + 
                      I*Sqrt[-3 - e^2 + p])), (4*e)/(-6 + 2*e + p)]]) - 
              2*Sqrt[-4*e^2 + (-2 + p)^2]*Pi*(-4*(-2 + 2*e + p)*EllipticF[
                  \[Chi]0, (4*e)/(-6 + 2*e + p)] - ((-1 + e)*(6 + 2*e - p)*
                  (-2 + 2*e + p)*Im[EllipticPi[(2*e*(-4 + p + (2*I)*
                        Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*(1 + e + 
                       I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/(-6 + 2*e + 
                      p)]])/Sqrt[-3 - e^2 + p] + (4*(1 + e)^2 - (-4 + p)^2)*
                 Re[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
                    ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), \[Chi]0, 
                   (4*e)/(-6 + 2*e + p)]])))/((4*e^2 - (-2 + p)^2)*
             (-2 + 2*e + p)*Pi)]*((1 + e)*(-6 + 2*e + p) - 
           2*e*(-4 + p)*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
              (4*e)/(-6 + 2*e + p)]^2))/((12 + 8*e - 4*e^2 - 8*p + p^2 - 
           16*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
              (4*e)/(-6 + 2*e + p)]^2)*Sqrt[(p*(-3 - e^2 + p)*
             (1 - ((3 + e^2)*(-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[
                        (4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]^2)^
                 2)/((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*
                  JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
                    (4*e)/(-6 + 2*e + p)]^2)^2 + (p*(-6 + 2*e + p - 
                  4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
                     (4*e)/(-6 + 2*e + p)]^2)^2)/((1 + e)*(-6 + 2*e + p) - 
                 2*e*(-4 + p)*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/
                     Pi, (4*e)/(-6 + 2*e + p)]^2)^2))/(-4*e^2 + 
             (-2 + p)^2)]) + (2*e*Abs[-4*e^2 + (-6 + p)^2]*
          Sqrt[(JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/
                (-6 + 2*e + p)]^2 - JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + 
                    p)])/Pi, (4*e)/(-6 + 2*e + p)]^4)/((3 + e^2 - p)*
             (-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + 
                       p)])/Pi, (4*e)/(-6 + 2*e + p)]^2)^2*(-12 - 8*e + 
              4*e^2 + 8*p - p^2 + 16*e*JacobiSN[(qr*EllipticK[(4*e)/
                     (-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]^2))]*
          Sin[(Sqrt[(-3 - e^2 + p)/(-6 + 2*e + p)]*
             ((2 + 2*e - p)*(-2 + 2*e + p)^2*Sqrt[(-6 + 2*e + p)/(-3 - e^2 + 
                  p)]*Pi*q\[Psi] + 2*Sqrt[-4*e^2 + (-2 + p)^2]*qr*(
                -4*(-2 + 2*e + p)*EllipticK[(4*e)/(-6 + 2*e + p)] - 
                ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                    (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                      (1 + e + I*Sqrt[-3 - e^2 + p])), (4*e)/(-6 + 2*e + 
                      p)]])/Sqrt[-3 - e^2 + p] + (4*(1 + e)^2 - (-4 + p)^2)*
                 Re[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
                    ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), 
                   (4*e)/(-6 + 2*e + p)]]) - 2*Sqrt[-4*e^2 + (-2 + p)^2]*Pi*(
                -4*(-2 + 2*e + p)*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)] - 
                ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                    (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                      (1 + e + I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/
                     (-6 + 2*e + p)]])/Sqrt[-3 - e^2 + p] + 
                (4*(1 + e)^2 - (-4 + p)^2)*Re[EllipticPi[(2*e*(-4 + p + 
                      (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*(1 + e + 
                      I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/(-6 + 2*e + 
                     p)]])))/((4*e^2 - (-2 + p)^2)*(-2 + 2*e + p)*Pi)])/
         Sqrt[(p*((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*JacobiSN[
                (qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + 
                  p)]^2)*((-3 - e^2 + p)^(-1) + (-6 + 2*e + p - 
                4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
                   (4*e)/(-6 + 2*e + p)]^2)^2/((1 + e)*(-6 + 2*e + p) - 
                2*e*(-4 + p)*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/
                    Pi, (4*e)/(-6 + 2*e + p)]^2)^2)*
            (-2 + (p*(-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[(4*e)/
                       (-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]^2))/
              ((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*JacobiSN[
                  (qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + 
                    p)]^2)))/(-6 + 2*e + p - 
            4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/
                (-6 + 2*e + p)]^2)]))/((1 + e)*(-6 + 2*e + p) - 
       2*e*(-4 + p)*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
          (4*e)/(-6 + 2*e + p)]^2), 
     (Sqrt[s^2 - sp^2]*Cos[(Sqrt[(-3 - e^2 + p)/(-6 + 2*e + p)]*
          ((2 + 2*e - p)*(-2 + 2*e + p)^2*Sqrt[(-6 + 2*e + p)/(-3 - e^2 + p)]*
            Pi*q\[Psi] + 2*Sqrt[-4*e^2 + (-2 + p)^2]*qr*
            (-4*(-2 + 2*e + p)*EllipticK[(4*e)/(-6 + 2*e + p)] - 
             ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                 (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                   (1 + e + I*Sqrt[-3 - e^2 + p])), (4*e)/(-6 + 2*e + p)]])/
              Sqrt[-3 - e^2 + p] + (4*(1 + e)^2 - (-4 + p)^2)*
              Re[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
                 ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), 
                (4*e)/(-6 + 2*e + p)]]) - 2*Sqrt[-4*e^2 + (-2 + p)^2]*Pi*
            (-4*(-2 + 2*e + p)*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)] - 
             ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                 (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                   (1 + e + I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/
                  (-6 + 2*e + p)]])/Sqrt[-3 - e^2 + p] + 
             (4*(1 + e)^2 - (-4 + p)^2)*Re[EllipticPi[(2*e*(-4 + p + 
                   (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*(1 + e + 
                   I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/(-6 + 2*e + p)]])))/
         ((4*e^2 - (-2 + p)^2)*(-2 + 2*e + p)*Pi)]*
       ((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*
          JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
            (4*e)/(-6 + 2*e + p)]^2)^2*Sqrt[(-3 - e^2 + p)^(-1) + 
         (-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/
                Pi, (4*e)/(-6 + 2*e + p)]^2)^2/((1 + e)*(-6 + 2*e + p) - 
            2*e*(-4 + p)*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
               (4*e)/(-6 + 2*e + p)]^2)^2])/
      (p*(-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/
             Pi, (4*e)/(-6 + 2*e + p)]^2)^2), 
     (Sqrt[s^2 - sp^2]*((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*
         JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)]^2)*(2*e*Abs[-4*e^2 + (-6 + p)^2]*
         Cos[(Sqrt[(-3 - e^2 + p)/(-6 + 2*e + p)]*
            ((2 + 2*e - p)*(-2 + 2*e + p)^2*Sqrt[(-6 + 2*e + p)/(-3 - e^2 + 
                 p)]*Pi*q\[Psi] + 2*Sqrt[-4*e^2 + (-2 + p)^2]*qr*
              (-4*(-2 + 2*e + p)*EllipticK[(4*e)/(-6 + 2*e + p)] - 
               ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                   (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                     (1 + e + I*Sqrt[-3 - e^2 + p])), (4*e)/(-6 + 2*e + p)]])/
                Sqrt[-3 - e^2 + p] + (4*(1 + e)^2 - (-4 + p)^2)*
                Re[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
                   ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), 
                  (4*e)/(-6 + 2*e + p)]]) - 2*Sqrt[-4*e^2 + (-2 + p)^2]*Pi*
              (-4*(-2 + 2*e + p)*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)] - 
               ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                   (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                     (1 + e + I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/
                    (-6 + 2*e + p)]])/Sqrt[-3 - e^2 + p] + (4*(1 + e)^2 - 
                 (-4 + p)^2)*Re[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[
                       -3 - e^2 + p]))/((-6 + 2*e + p)*(1 + e + 
                     I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/(-6 + 2*e + 
                    p)]])))/((4*e^2 - (-2 + p)^2)*(-2 + 2*e + p)*Pi)]*
         Sqrt[(p*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
              (4*e)/(-6 + 2*e + p)]^2*(-1 + JacobiSN[(qr*EllipticK[
                  (4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]^2))/
           ((3 + e^2 - p)*(6 - 2*e - p + 4*e*JacobiSN[(qr*EllipticK[
                   (4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]^2)*
            (-((-6 + 2*e + p)^2*(-2 + 2*e + p)) - 4*e*(-60 - 4*e^2 + 4*e^3 - 
               e*(-2 + p)^2 + 28*p - 3*p^2)*JacobiSN[(qr*EllipticK[
                   (4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]^2 + 
             4*e^2*(4*e^2 - (-2 + p)^2)*JacobiSN[(qr*EllipticK[(4*e)/
                    (-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]^4))] + 
        (Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(-6 + 2*e + p - 
            4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/
                (-6 + 2*e + p)]^2)^2*Sin[(Sqrt[(-3 - e^2 + p)/(-6 + 2*e + p)]*
             ((2 + 2*e - p)*(-2 + 2*e + p)^2*Sqrt[(-6 + 2*e + p)/(-3 - e^2 + 
                  p)]*Pi*q\[Psi] + 2*Sqrt[-4*e^2 + (-2 + p)^2]*qr*(
                -4*(-2 + 2*e + p)*EllipticK[(4*e)/(-6 + 2*e + p)] - 
                ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                    (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                      (1 + e + I*Sqrt[-3 - e^2 + p])), (4*e)/(-6 + 2*e + 
                      p)]])/Sqrt[-3 - e^2 + p] + (4*(1 + e)^2 - (-4 + p)^2)*
                 Re[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
                    ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), 
                   (4*e)/(-6 + 2*e + p)]]) - 2*Sqrt[-4*e^2 + (-2 + p)^2]*Pi*(
                -4*(-2 + 2*e + p)*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)] - 
                ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                    (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                      (1 + e + I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/
                     (-6 + 2*e + p)]])/Sqrt[-3 - e^2 + p] + 
                (4*(1 + e)^2 - (-4 + p)^2)*Re[EllipticPi[(2*e*(-4 + p + 
                      (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*(1 + e + 
                      I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/(-6 + 2*e + 
                     p)]])))/((4*e^2 - (-2 + p)^2)*(-2 + 2*e + p)*Pi)])/
         (((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*
             JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/
                (-6 + 2*e + p)]^2)^2*
          Sqrt[1 - ((3 + e^2)*(-6 + 2*e + p - 4*e*JacobiSN[
                   (qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + 
                     p)]^2)^2)/((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*
                JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
                  (4*e)/(-6 + 2*e + p)]^2)^2 + 
            (p*(-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + 
                        p)])/Pi, (4*e)/(-6 + 2*e + p)]^2)^2)/
             ((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*JacobiSN[
                  (qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + 
                    p)]^2)^2])))/(p*(-6 + 2*e + p - 
        4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)]^2)), 
     (sp*((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*
         JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)]^2))/(p*(-6 + 2*e + p - 
        4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)]^2))}
 
\[Delta]t = -1/2*((-1 + e)*(3 + e)*(-4*e^2 + (-2 + p)^2)*p*qr*
        (-((p*(-6 + 2*e + p)*EllipticE[(4*e)/(-6 + 2*e + p)])/
           ((-1 + e)*(1 + e)*(-4 + p))) + (p*(36 + 2*e^2*(-2 + p) - 14*p + 
            p^2)*EllipticK[(4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^2) - 
         (2*(6 + 2*e - p)*(8 + p - p^2 + e^2*(-8 + 3*p))*
           EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), 
            (4*e)/(-6 + 2*e + p)])/((-1 + e)*(1 + e)^2*(-4 + p)^2) - 
         (2*(6 + 2*e - p)*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
            (4*e)/(-6 + 2*e + p)])/(2 + 2*e - p)))/
       ((-((p*(-6 + 2*e + p))/(3 + e^2 - p)))^(3/2)*(-3 - e^2 + p)^(5/2)*
        Pi) + ((-1 + e^2)^2*(-3 - e^2 + p)*qr*
       (-((p*(-6 + 2*e + p)*EllipticE[(4*e)/(-6 + 2*e + p)])/
          ((-1 + e)*(1 + e)*(-4 + p))) + 
        (p*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*EllipticK[
           (4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^2) - 
        (2*(6 + 2*e - p)*(8 + p - p^2 + e^2*(-8 + 3*p))*
          EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), 
           (4*e)/(-6 + 2*e + p)])/((-1 + e)*(1 + e)^2*(-4 + p)^2) - 
        (2*(6 + 2*e - p)*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
           (4*e)/(-6 + 2*e + p)])/(2 + 2*e - p)))/(2*(3 + e^2 - p)^2*
       Sqrt[p*(-6 + 2*e + p)]*Pi) - 
     (Sqrt[(-4*e^2 + (-2 + p)^2)/(p^2*(-6 + 2*e + p))]*qr*
       ((4*p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*EllipticE[(4*e)/(-6 + 2*e + p)])/
         ((1 + e)*(-4 + p)^2) + (p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*
          (EllipticE[(4*e)/(-6 + 2*e + p)] - EllipticK[
            (4*e)/(-6 + 2*e + p)]))/((-1 + e)*(1 + e)*(-4 + p)) + 
        (p^2*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*
          EllipticK[(4*e)/(-6 + 2*e + p)])/(2*(-1 + e^2)*(-4 + p)^3) - 
        (p^(3/2)*(128 - 144*p - 32*p^2 + 64*p^3 - 17*p^4 + p^5 - 
           16*e^4*(-8 + 5*p) + e^2*(-256 - 32*p + 32*p^2 + p^4))*
          EllipticK[(4*e)/(-6 + 2*e + p)])/(2*(-1 + e^2)*
          Sqrt[-4*e^2 + (-2 + p)^2]*(-4 + p)^3) + 
        (p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*
          ((-6 + 2*e + p)*EllipticE[(4*e)/(-6 + 2*e + p)] + 
           (6 + 2*e - p)*EllipticK[(4*e)/(-6 + 2*e + p)]))/
         ((-1 + e)*(1 + e)*(6 + 2*e - p)*(-4 + p)^2*(-6 + 2*e + p)) - 
        (4*p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(-6 - 2*e + p)*
          EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), 
           (4*e)/(-6 + 2*e + p)])/((1 + e)*(-4 + p)^3) - 
        (8*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(-8 + e^2*(8 - 3*p) - p + p^2)*
          EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), 
           (4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^3) - 
        (4*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(-8 + e^2*(8 - 3*p) - p + p^2)*
          ((1 + e)*(-6 + 2*e + p)*EllipticK[(4*e)/(-6 + 2*e + p)] - 
           (-6 + 2*e^2 + 3*e*(-4 + p) + p)*EllipticPi[(2*e*(-4 + p))/
              ((1 + e)*(-6 + 2*e + p)), (4*e)/(-6 + 2*e + p)]))/
         ((-1 + e)*(1 + e)^2*(-4 + p)^3*(-6 + 2*e + p)) - 
        (8*(1 + e)*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*EllipticPi[
           (16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), (4*e)/(-6 + 2*e + p)])/
         ((-4 + p)*(-2 - 2*e + p)) + (p^(5/2)*(-6 - 2*e + p)*
          EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
           (4*e)/(-6 + 2*e + p)])/(2*Sqrt[-4*e^2 + (-2 + p)^2]*
          (-2 - 2*e + p)) - (p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(-6 - 2*e + p)*
          EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
           (4*e)/(-6 + 2*e + p)])/((-4 + p)*(-2 - 2*e + p)) + 
        (Sqrt[(-4*e^2 + (-2 + p)^2)*p]*((12 + 8*e - 4*e^2 - 8*p + p^2)*
            EllipticK[(4*e)/(-6 + 2*e + p)] + (-12 - 24*e + 4*e^2 + 8*p - 
             p^2)*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
             (4*e)/(-6 + 2*e + p)]))/(2*(2 + 2*e - p)*(-6 + 2*e + p))))/
      (2*Pi) + ((-1 + e)*(3 + e)*(-4*e^2 + (-2 + p)^2)*p*
       (-((p*(-6 + 2*e + p)*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)])/
          ((-1 + e)*(1 + e)*(-4 + p))) + 
        (p*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*EllipticF[\[Chi]0, 
           (4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^2) - 
        (2*(6 + 2*e - p)*(8 + p - p^2 + e^2*(-8 + 3*p))*
          EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), \[Chi]0, 
           (4*e)/(-6 + 2*e + p)])/((-1 + e)*(1 + e)^2*(-4 + p)^2) - 
        (2*(6 + 2*e - p)*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
           \[Chi]0, (4*e)/(-6 + 2*e + p)])/(2 + 2*e - p) + 
        (e*p*Sqrt[(-6 + 2*e + p)*(-6 + p + 2*e*Cos[2*\[Chi]0])]*
          Sin[2*\[Chi]0])/((-1 + e^2)*(-6 + 2*e^2 + p + 
           e*(-4 + p)*Cos[2*\[Chi]0]))))/
      (2*(-((p*(-6 + 2*e + p))/(3 + e^2 - p)))^(3/2)*(-3 - e^2 + p)^(5/2)) - 
     ((-1 + e^2)^2*(-3 - e^2 + p)*
       (-((p*(-6 + 2*e + p)*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)])/
          ((-1 + e)*(1 + e)*(-4 + p))) + 
        (p*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*EllipticF[\[Chi]0, 
           (4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^2) - 
        (2*(6 + 2*e - p)*(8 + p - p^2 + e^2*(-8 + 3*p))*
          EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), \[Chi]0, 
           (4*e)/(-6 + 2*e + p)])/((-1 + e)*(1 + e)^2*(-4 + p)^2) - 
        (2*(6 + 2*e - p)*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), 
           \[Chi]0, (4*e)/(-6 + 2*e + p)])/(2 + 2*e - p) + 
        (e*p*Sqrt[(-6 + 2*e + p)*(-6 + p + 2*e*Cos[2*\[Chi]0])]*
          Sin[2*\[Chi]0])/((-1 + e^2)*(-6 + 2*e^2 + p + 
           e*(-4 + p)*Cos[2*\[Chi]0]))))/(2*(3 + e^2 - p)^2*
       Sqrt[p*(-6 + 2*e + p)]) + 
     Sqrt[(-4*e^2 + (-2 + p)^2)/(p^2*(-6 + 2*e + p))]*
      ((2*p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*EllipticE[\[Chi]0, 
          (4*e)/(-6 + 2*e + p)])/((1 + e)*(-4 + p)^2) + 
       (p^2*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*
         EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)])/(4*(-1 + e^2)*
         (-4 + p)^3) - (p^(3/2)*(128 - 144*p - 32*p^2 + 64*p^3 - 17*p^4 + 
          p^5 - 16*e^4*(-8 + 5*p) + e^2*(-256 - 32*p + 32*p^2 + p^4))*
         EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)])/(4*(-1 + e^2)*
         Sqrt[-4*e^2 + (-2 + p)^2]*(-4 + p)^3) + 
       (p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(36 + 2*e^2*(-2 + p) - 14*p + p^2)*
         (qr*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*Cos[2*\[Chi]0]]*
           EllipticE[(4*e)/(-6 + 2*e + p)] + 
          Pi*(-6*Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/(-6 + 2*e + p)] + 
            2*e*Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/(-6 + 2*e + p)] + 
            p*Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/(-6 + 2*e + p)] - 
            Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*Cos[2*\[Chi]0]])*
           EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)] + (6 + 2*e - p)*Pi*
           Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/(-6 + 2*e + p)]*
           EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)]))/(2*(-1 + e)*(1 + e)*
         (6 + 2*e - p)*(-4 + p)^2*Pi*Sqrt[(-6 + 2*e + p)*
           (-6 + p + 2*e*Cos[2*\[Chi]0])]) - 
       (2*p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(-6 - 2*e + p)*
         EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), \[Chi]0, 
          (4*e)/(-6 + 2*e + p)])/((1 + e)*(-4 + p)^3) - 
       (4*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(-8 + e^2*(8 - 3*p) - p + p^2)*
         EllipticPi[(2*e*(-4 + p))/((1 + e)*(-6 + 2*e + p)), \[Chi]0, 
          (4*e)/(-6 + 2*e + p)])/((-1 + e^2)*(-4 + p)^3) - 
       (4*(1 + e)*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*EllipticPi[
          (16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), \[Chi]0, 
          (4*e)/(-6 + 2*e + p)])/((-4 + p)*(-2 - 2*e + p)) + 
       (p^(5/2)*(-6 - 2*e + p)*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + 
            p^2), \[Chi]0, (4*e)/(-6 + 2*e + p)])/
        (4*Sqrt[-4*e^2 + (-2 + p)^2]*(-2 - 2*e + p)) - 
       (p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(-6 - 2*e + p)*
         EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + p^2), \[Chi]0, 
          (4*e)/(-6 + 2*e + p)])/(2*(-4 + p)*(-2 - 2*e + p)) + 
       (e*p^(3/2)*Cos[2*\[Chi]0]*Sqrt[((-4*e^2 + (-2 + p)^2)*
            (-6 + p + 2*e*Cos[2*\[Chi]0]))/(-6 + 2*e + p)]*
         (qr*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*Cos[2*\[Chi]0]]*
           EllipticE[(4*e)/(-6 + 2*e + p)] - 
          Pi*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*Cos[2*\[Chi]0]]*
           EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)] + 2*e*Pi*Sin[2*\[Chi]0]))/
        ((-1 + e^2)*(6 + 2*e - p)*Pi*(-6 + 2*e^2 + p + 
          e*(-4 + p)*Cos[2*\[Chi]0])) - (p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*
         (qr*Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/(-6 + 2*e + p)]*
           Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*Cos[2*\[Chi]0]]*
           EllipticE[(4*e)/(-6 + 2*e + p)] + 
          Pi*((-6 - 2*e + p - Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/(-6 + 2*e + 
                  p)]*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 
                 8*e^2*Cos[2*\[Chi]0]])*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + 
                p)] + (6 + 2*e - p)*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + 
                p)] + 2*e*Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/(-6 + 2*e + p)]*
             Sin[2*\[Chi]0])))/(2*(-1 + e)*(1 + e)*(6 + 2*e - p)*(-4 + p)*
         Pi) - ((1 + e)*((1 + e)^(-1) - 2/(-4 + p))*(-4 + p)*
         Sqrt[(-4*e^2 + (-2 + p)^2)*p]*
         ((4*(2 + 2*e - p)*(-(Sqrt[(-6 + 2*e + p)*(-6 + p + 
                  2*e*Cos[2*\[Chi]0])]*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + 
                  p)]) - (6 + 2*e - p)*Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/
                (-6 + 2*e + p)]*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 8*p + 
                 p^2), \[Chi]0, (4*e)/(-6 + 2*e + p)] + 2*e*Sin[2*\[Chi]0]))/
           ((-6 - 2*e + p)^2*Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/
              (-6 + 2*e + p)]) + (4*(2 + 2*e - p)*(-6 + 2*e + p)^(3/2)*
            (qr*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*
                 Cos[2*\[Chi]0]]*EllipticE[(4*e)/(-6 + 2*e + p)] - 
             Pi*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*
                 Cos[2*\[Chi]0]]*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)] + 
             2*e*Pi*Sin[2*\[Chi]0]))/((6 + 2*e - p)*Pi*
            Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]]*(12 - 4*e^2 - 8*p + p^2 + 
             8*e*Cos[2*\[Chi]0])) - (4*(2 + 2*e - p)*(-6 + 2*e + p)*
             (-12 + 4*e^2 + 8*p - p^2 - 8*e*Cos[2*\[Chi]0])*
             EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)] - (2 + 2*e - p)*
             (6 + 2*e - p)*(-6 + 2*e + p)*(-12 + 4*e^2 + 8*p - p^2 - 
              8*e*Cos[2*\[Chi]0])*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)] + 
            (64*e - (-2 - 2*e + p)^2*(-6 + 2*e + p))*(12 - 4*e^2 - 8*p + 
              p^2 + 8*e*Cos[2*\[Chi]0])*EllipticPi[(16*e)/(12 + 8*e - 4*e^2 - 
                8*p + p^2), \[Chi]0, (4*e)/(-6 + 2*e + p)] + 
            32*e*(2 + 2*e - p)*Sqrt[(-6 + 2*e + p)*(-6 + p + 
                2*e*Cos[2*\[Chi]0])]*Sin[2*\[Chi]0])/((-6 - 2*e + p)^2*
            (12 - 4*e^2 - 8*p + p^2 + 8*e*Cos[2*\[Chi]0]))))/
        (4*(2 + 2*e - p)*(-6 + 2*e + p)) - (((1 + e)^(-1) - 2/(-4 + p))*
         Sqrt[(-4*e^2 + (-2 + p)^2)*p]*(4 + ((1 - e)^(-1) + (1 + e)^(-1) + 
            2/(-4 + p))*p)*((2*(1 + e)*(-(Sqrt[(-6 + 2*e + p)*(-6 + p + 
                  2*e*Cos[2*\[Chi]0])]*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + 
                  p)]) - (6 + 2*e - p)*Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/
                (-6 + 2*e + p)]*EllipticPi[(2*e*(-4 + p))/((1 + e)*
                 (-6 + 2*e + p)), \[Chi]0, (4*e)/(-6 + 2*e + p)] + 
             2*e*Sin[2*\[Chi]0]))/Sqrt[(-6 + p + 2*e*Cos[2*\[Chi]0])/
             (-6 + 2*e + p)] - ((1 + e)*(6 + 2*e - p)*(-6 + 2*e + p)^(3/2)*
            (qr*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*
                 Cos[2*\[Chi]0]]*EllipticE[(4*e)/(-6 + 2*e + p)] - 
             Pi*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*
                 Cos[2*\[Chi]0]]*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)] + 
             2*e*Pi*Sin[2*\[Chi]0]))/(Pi*Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]]*
            (-6 + 2*e^2 + p + e*(-4 + p)*Cos[2*\[Chi]0])) + 
          (2*(1 - e)*(-((1 + e)*(-4 + p)*(-6 + 2*e + p)*(-6 + 2*e^2 + p + 
                e*(-4 + p)*Cos[2*\[Chi]0])*EllipticE[\[Chi]0, 
                (4*e)/(-6 + 2*e + p)]) - (1 + e)*(6 + 2*e - p)*(-6 + 2*e + p)*
              (-6 + 2*e^2 + p + e*(-4 + p)*Cos[2*\[Chi]0])*EllipticF[\[Chi]0, 
               (4*e)/(-6 + 2*e + p)] + 2*(-6 + 2*e^3 + e^2*(-2 + p) + p - e*
                (26 - 10*p + p^2))*(-6 + 2*e^2 + p + e*(-4 + p)*
                Cos[2*\[Chi]0])*EllipticPi[(2*e*(-4 + p))/((1 + e)*
                 (-6 + 2*e + p)), \[Chi]0, (4*e)/(-6 + 2*e + p)] + 
             e*(1 + e)*(-4 + p)^2*Sqrt[(-6 + 2*e + p)*(-6 + p + 
                 2*e*Cos[2*\[Chi]0])]*Sin[2*\[Chi]0]))/((-1 + e)*(-4 + p)*
            (-6 + 2*e^2 + p + e*(-4 + p)*Cos[2*\[Chi]0]))))/
        (2*(-6 - 2*e + p)^2*(-6 + 2*e + p)) + 
       (e*p*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*Sin[2*\[Chi]0]*
         (-(e*(-4 + p)*qr*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 
              8*e^2*Cos[2*\[Chi]0]]*EllipticE[(4*e)/(-6 + 2*e + p)]*
            Sin[2*\[Chi]0]) + Pi*(-72 + 24*e + 20*e^2 + 24*p - 10*e*p - 
            3*e^2*p - 2*p^2 + e*p^2 - e*(-48 + 8*e^2 - 2*e*(-4 + p) + 14*p - 
              p^2)*Cos[2*\[Chi]0] - 4*e^2*Cos[4*\[Chi]0] + 
            e^2*p*Cos[4*\[Chi]0] + e*(-4 + p)*Sqrt[-4*e^2 + 4*e*(-6 + p) + 
               (-6 + p)^2 + 8*e^2*Cos[2*\[Chi]0]]*EllipticE[\[Chi]0, 
              (4*e)/(-6 + 2*e + p)]*Sin[2*\[Chi]0])))/((-1 + e)*(1 + e)*
         (6 + 2*e - p)*(-4 + p)*Pi*Sqrt[(-6 + 2*e + p)*
           (-6 + p + 2*e*Cos[2*\[Chi]0])]*(-6 + 2*e^2 + p + 
          e*(-4 + p)*Cos[2*\[Chi]0])) - 
       (e*p^(3/2)*Sqrt[((-4*e^2 + (-2 + p)^2)*(-6 + p + 2*e*Cos[2*\[Chi]0]))/
           (-6 + 2*e + p)]*Sin[2*\[Chi]0]*
         (-(e*(-4 + p)^2*qr*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 
              8*e^2*Cos[2*\[Chi]0]]*EllipticE[(4*e)/(-6 + 2*e + p)]*
            Sin[2*\[Chi]0]) + Pi*(-72 + 64*e^2 - 8*e^4 + 24*p - 16*e^2*p - 
            2*p^2 + e^2*p^2 + e^2*(-4 + p)^2*Cos[4*\[Chi]0] + 
            e*(-4 + p)^2*Sqrt[-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 8*e^2*
                Cos[2*\[Chi]0]]*EllipticE[\[Chi]0, (4*e)/(-6 + 2*e + p)]*
             Sin[2*\[Chi]0])))/((-1 + e^2)*(6 + 2*e - p)*(-4 + p)*Pi*
         (-6 + 2*e^2 + p + e*(-4 + p)*Cos[2*\[Chi]0])^2))
 
\[Delta]r = (-2*e*Sqrt[(-4*e^2 + (-2 + p)^2)*p]*
      JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]*
      (JacobiCN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
         (4*e)/(-6 + 2*e + p)]*JacobiDN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/
          Pi, (4*e)/(-6 + 2*e + p)]*((-6 + 2*e + p)*qr*
          EllipticE[(4*e)/(-6 + 2*e + p)] - (-6 + 2*e + p)*Pi*
          JacobiEpsilon[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)] + 4*e*Pi*JacobiCD[
           (qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]*
          JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)]) + 4*e*Pi*JacobiSN[
         (qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]*
        (-1 + JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)]^2)))/
     (Pi*((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*
         JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
           (4*e)/(-6 + 2*e + p)]^2)^2)
 
\[Delta]\[CurlyPhi] = (I*Sqrt[4*e^2 - (-2 + p)^2]*
      (p*qr*(-6*Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]] + 
         2*e*Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]] + 
         p*Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]] - 
         Sqrt[(-6 + 2*e + p)*(-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 
            8*e^2*Cos[2*\[Chi]0])])*EllipticE[(4*e)/(-6 + 2*e + p)] - 
       p*Pi*(-6*Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]] + 
         2*e*Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]] + 
         p*Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]] - 
         Sqrt[(-6 + 2*e + p)*(-4*e^2 + 4*e*(-6 + p) + (-6 + p)^2 + 
            8*e^2*Cos[2*\[Chi]0])])*EllipticE[\[Chi]0, 
         (4*e)/(-6 + 2*e + p)] + (-4*e^2 + (-6 + p)^2)*
        Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]]*
        (Pi*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)] - 
         qr*EllipticK[(4*e)/(-6 + 2*e + p)])))/((6 + 2*e - p)*p*
      (-6 + 2*e + p)^(3/2)*Pi*Sqrt[-6 + p + 2*e*Cos[2*\[Chi]0]])
 
\[Delta]\[CurlyTheta] = (((1 + e)*(-6 + 2*e + p) - 
       2*e*(-4 + p)*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
          (4*e)/(-6 + 2*e + p)]^2)*Sqrt[(s - sp)*(s + sp)*
        (1 - ((3 + e^2)*(-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[
                   (4*e)/(-6 + 2*e + p)])/Pi, (4*e)/(-6 + 2*e + p)]^2)^2)/
          ((1 + e)*(-6 + 2*e + p) - 2*e*(-4 + p)*
             JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, (4*e)/
                (-6 + 2*e + p)]^2)^2 + 
         (p*(-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/
                 Pi, (4*e)/(-6 + 2*e + p)]^2)^2)/((1 + e)*(-6 + 2*e + p) - 
            2*e*(-4 + p)*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
               (4*e)/(-6 + 2*e + p)]^2)^2)]*
      Sin[(Sqrt[(-3 - e^2 + p)/(-6 + 2*e + p)]*
         ((2 + 2*e - p)*(-2 + 2*e + p)^2*Sqrt[(-6 + 2*e + p)/(-3 - e^2 + p)]*
           Pi*q\[Psi] + 2*Sqrt[-4*e^2 + (-2 + p)^2]*qr*
           (-4*(-2 + 2*e + p)*EllipticK[(4*e)/(-6 + 2*e + p)] - 
            ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                  (1 + e + I*Sqrt[-3 - e^2 + p])), (4*e)/(-6 + 2*e + p)]])/
             Sqrt[-3 - e^2 + p] + (4*(1 + e)^2 - (-4 + p)^2)*
             Re[EllipticPi[(2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/
                ((-6 + 2*e + p)*(1 + e + I*Sqrt[-3 - e^2 + p])), (4*e)/
                (-6 + 2*e + p)]]) - 2*Sqrt[-4*e^2 + (-2 + p)^2]*Pi*
           (-4*(-2 + 2*e + p)*EllipticF[\[Chi]0, (4*e)/(-6 + 2*e + p)] - 
            ((-1 + e)*(6 + 2*e - p)*(-2 + 2*e + p)*Im[EllipticPi[
                (2*e*(-4 + p + (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*
                  (1 + e + I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/
                 (-6 + 2*e + p)]])/Sqrt[-3 - e^2 + p] + 
            (4*(1 + e)^2 - (-4 + p)^2)*Re[EllipticPi[(2*e*(-4 + p + 
                  (2*I)*Sqrt[-3 - e^2 + p]))/((-6 + 2*e + p)*(1 + e + 
                  I*Sqrt[-3 - e^2 + p])), \[Chi]0, (4*e)/(-6 + 2*e + p)]])))/
        ((4*e^2 - (-2 + p)^2)*(-2 + 2*e + p)*Pi)])/
     (p*(-6 + 2*e + p - 4*e*JacobiSN[(qr*EllipticK[(4*e)/(-6 + 2*e + p)])/Pi, 
          (4*e)/(-6 + 2*e + p)]^2))
