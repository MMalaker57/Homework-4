//
//  Schrodinger RK4.swift
//  Test Plot Threaded
//
//  Created by Matthew Malaker on 2/25/22.
//

import Foundation
import SwiftUI
import CorePlot


class Schrodinger: NSObject, ObservableObject {
 
    var potential = Potential()
//    var potentialV: (onedArray: [Double], xArray: [Double], yArray: [Double]) = ([],[],[])
    
    
    var hbar2over2m = 0.0
    
    //k and l are the same type of variable but are used in the separate RK4 solutions
    //we need 2 solutions because we have a second order equation and need to solve for psi' and psi, not just psi
    //We do this because we only have a relation for psi' from psi'', which itself is a variation of psi, but psi' is unknown
    //past the initial condition.
//    var k1: Double = 0.0
//    var k2: Double = 0.0
//    var k3: Double = 0.0
//    var k4: Double = 0.0
//    var l1: Double = 0.0
//    var l2: Double = 0.0
//    var l3: Double = 0.0
//    var l4: Double = 0.0
    
    
//      2        2
//- hbar  partial psi
//-------- ----------- + VPsi = EPsi
//   2m              2
//          partial x

    

    
    
//      2        2
//- hbar  partial psi
//-------- ----------- + (V-E)Psi = 0
//   2m              2
//          partial x

    
    
// We need to solve for E, Psi, and the functional
    
//
//Psi'' = (2m(V-E)/hbar^2)*Psi
//
//
//
//
    
    override init() {
        
        
        //Must call super init before initializing plot
        super.init()
        
        let hbar = 6.582119569e-16
        let hbar2 = pow((6.582119569e-16),2.0)
        
        let massE = 510998.946/(pow(299792458.0e10,2.0))
        
        hbar2over2m = hbar2/(2.0*massE)
//      print(hbar2over2m)
    }

    func setV(choice: String, xMin: Double, xMax: Double, stepSize: Double){
        potential.getPotential(potentialToGet: choice, xMin: xMin, xMax: xMax, stepSize: stepSize)
        
    }
    
    func calculatePsiForE(xMin: Double, xMax: Double, stepSize: Double, potentialStr: String, E: Double)->[Double]{
        
        var potentialString = ""
        var psi: [Double] = []
        var psiforE: [(E: Double, psi: [Double])] = []
        var psiPrime: [Double] = []
        var psiDoublePrime: [Double] = []
        
        
        var psiWithE: [Double] = []
        psi.removeAll()
        psiPrime.removeAll()
        psi.append(0.0)
        psiPrime.append(5.0)
        
//        potential.getPotential(potentialToGet: potentialStr, xMin: xMin, xMax: xMax, stepSize: stepSize)
        var potentialV = potential.potential
        
//        var eStep = 0.1
//        var maxE = 30.0
        var maeta = Int((xMax-xMin)/stepSize)
        
        
        
        //Solving a second order equation requires we solve two first orders. We have that

        var currentPsi = 0.0
        var nextPsi = 0.0
        var currentPsiPrime = psiPrime[0]
        var nextPsiPrime = 0.0
        var currentPsiDoublePrime = 0.0
        
    
            for i in stride(from: 0, to: maeta, by: 1){
               
//                print("V = \(potentialV.yArray[i])")
                currentPsiDoublePrime = ((potentialV.yArray[i]-E)*currentPsi)/hbar2over2m
                
                let k1 = currentPsiPrime
                let l1 = currentPsiDoublePrime
                
                
                let VAverage = ((potentialV.yArray[i+1]-potentialV.yArray[i])/2.0)
                
                
                let k2 = currentPsiPrime + (stepSize * (k1 / 2.0))
                let l2 = ((VAverage-E)/hbar2over2m) * (currentPsi + (stepSize*l1)/2.0)

                let k3 = currentPsiPrime + (stepSize * (k2 / 2.0))
                let l3 = ((VAverage-E)/hbar2over2m) * (currentPsi + (stepSize*l2)/2.0)

                let k4 = currentPsiPrime + (stepSize * k3)
                let l4 = ((potentialV.yArray[i+1]-E)/hbar2over2m) * (currentPsi + (stepSize*l3))
                
        
                nextPsi = currentPsi + (stepSize/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
                nextPsiPrime = currentPsiPrime + (stepSize/6.0)*(l1 + 2.0*l2 + 2.0*l3 + l4)
                
//                nextPsi = currentPsi + stepSize*k1
//                nextPsiPrime = currentPsiPrime + stepSize*l1
     
                
//                print("Psi = \(nextPsi)")
//                print("PsiPrime = \(nextPsiPrime)")
                psi.append(nextPsi)
                psiPrime.append(nextPsiPrime)
                currentPsi = nextPsi
                currentPsiPrime = nextPsiPrime
                
            }
        
            nextPsi = 0.0
            nextPsiPrime = 0.0
            currentPsi = 0.0
            currentPsiPrime = 0.0
            potentialString = ""
            psiforE = []
            psiPrime = []
            psiDoublePrime = []
        
            psiWithE = psi
        
        return psiWithE
        
        
        
    }

    func separateGoodE(minX: Double, maxX: Double, stepSize: Double, potentialStr: String, data: [(E: Double, psi: [Double])])->[(E: Double, psi: [Double])] {
        var psiforESep: [(E: Double, psi: [Double])] = []
        var EList: [Double] = []
        var psiAtA: [Double] = []
        var psiList: [[Double]] = []
        let maeta = Int((maxX - minX) / stepSize)
        var pointsBeforeCross: [Int] = []
        var numberofZeroes: Int = 0
        
        let twoStep = 2.0 * stepSize //Just to not double stepSize many times. It's best to do it once
        
        
        for i in data{
            var temp = i.psi.last ?? 0.0
            psiAtA.append(temp) //If you put he definition of temp into this line, you get an error.
            psiList.append(i.psi)
            EList.append(i.E)
        }
        
        //We now have a function psi(A,E). We need to calculate the zeroes
        
        //To do this, we need the derivative of the function at each point we use. We can either calculate it at each point
        //or calculate it as we go. The latter is more efficient
        
        
        //We first should find how many zeroes there are
        //ironically, this is not a terrible way to find the zeroes

        for i in stride(from: 1, to: psiAtA.count-1, by: 1){
            if(psiAtA[i]*psiAtA[i+1] < 0){
                numberofZeroes += 1
                pointsBeforeCross.append(i)
            }
            
        }
        
        //Calculate zeroes
        for root in stride(from: 0, to: numberofZeroes, by: 1){
        
            var eta = pointsBeforeCross[root]
            var calculatedPsi: [Double] = []
            var E = Double(EList[pointsBeforeCross[root]])
            var functionalPrime = (psiAtA[eta+1] - psiAtA[eta])/((EList[eta+1]-EList[eta])) // initial
            var currentValue = psiAtA[eta]
            for i in stride(from: 1, to: 20, by: 1){
                   
                //Break if the derivative is too flat and continuing would cause problem
//                if(functionalPrime < 0.01){
//                    psiforESep.append((E: EList[eta], psi: psiList[eta]))
//                    break
//                }
                
                
                print("E= \(E) and currentValue = \(currentValue) and functionalPrime = \(functionalPrime)")
                let oldE = E
                E = oldE - (currentValue/functionalPrime)
            
                
                let oldFunctionalValue = currentValue
                calculatedPsi = calculatePsiForE(xMin: minX, xMax: maxX, stepSize: stepSize, potentialStr: potentialStr, E: E)
                currentValue = calculatedPsi.last!
                
                
                //derivative
                functionalPrime = (currentValue - oldFunctionalValue)/((E-oldE))
                
                if(abs(currentValue) < 0.0001){
                    psiforESep.append((E: E, psi: calculatedPsi))
                    break
                }
            
            }
            
        }
    
        return psiforESep
    }
    
    func normalize(data: (E: Double, psi: [Double]), stepSize: Double) -> (E: Double, psi: [Double]){
        var psiarray = data.psi
        var normarray: [Double] = []
        var normed: (E: Double, psi: [Double]) = (E: data.E, [])
        var integral: Double = 0.0
        
        //We need to integrate the two wavefunctions
        
        for i in psiarray{
            integral += pow(i,2.0)*stepSize
        }
        for i in psiarray{
            normarray.append(i*pow(integral,-0.5))
        }
        
        normed = (data.E, normarray)
        return normed
    }
    
}
