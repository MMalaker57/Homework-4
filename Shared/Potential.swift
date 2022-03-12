//
//  RK4 Schrodinger.swift
//  Test Plot Threaded
//
//  Created by Matthew Malaker on 2/25/22.
//

import Foundation
import SwiftUI
import CorePlot



class Potential: NSObject, ObservableObject{
    
    var potential: (onedArray: [Double], xArray: [Double], yArray: [Double]) = (onedArray: [],xArray: [] ,yArray: [])
    var xOffset = 0.0
    var hbar2over2m: Double = 0.0
    var hbar2overm: Double = 0.0
    var contentArray:[(xPoint: Double, yPoint: Double)] = []
    override init() {
        
        
        //Must call super init before initializing plot
        super.init()
        
        let hbar = 6.582119569e-16
        let hbar2 = pow((6.582119569e-16),2.0)
        
        let massE = 510998.946/(pow(299792458.0e10,2.0))
        
        hbar2over2m = hbar2/(2.0*massE)
        hbar2overm = hbar2/massE

    }
    
    
    
    func getPotential(potentialToGet: String, xMin: Double, xMax: Double, stepSize: Double){
        
        print(" Potential To Get is: \(potentialToGet)")
        potential.onedArray = []
        potential.xArray = []
        potential.yArray = []
        var count = 0
        var dataPoint: (xPoint: Double, yPoint: Double)
        contentArray.removeAll()

        switch potentialToGet{
            
        case "Square Well":
            print("Getting Square Well")
            //We set the bounds of the potential as "infinite" and the rest as zero
            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
            //Potential is in index space, not x space
            let maxIndex = Int((xMax-stepSize)/stepSize)
            for i in stride(from: 1, through: maxIndex, by: 1){
                print(i)
                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append(0.0)
            }
            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
            
            
        case "Linear Well":
            print("Getting Linear Well")
            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
            let maxIndex = Int((xMax-stepSize)/stepSize)
            for i in stride(from: 1, through: maxIndex, by: 1) {

                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append((Double(i)*stepSize+xMin-xMin)*4.0*1.3)
                //potential.yArray.append((i-xMin)*0.25)

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)

                
            }
            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
        
            
        case "Parabolic Well":
            print("Getting Parabolic Well")
            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
            let maxIndex = Int((xMax-stepSize)/stepSize)
        
            for i in stride(from: 1, through: maxIndex, by: 1) {

                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append((pow(((Double(i)*stepSize+xMin)-(xMax+xMin)/2.0), 2.0)/1.0))
                print(((pow(((Double(i)*stepSize+xMin)-(xMax+xMin)/2.0), 2.0)/1.0)))
                
                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)
                
            }
        
            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
        
        case "Square + Linear Well":
            print("Getting Square + Linear Well")
                    
            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
            let maxIndex = Int((xMax-stepSize)/stepSize)
            for i in stride(from: 1, to: Int((maxIndex)/2), by: 1) {

                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append(0.0)

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)

            }

            for i in stride(from: Int((maxIndex)/2), through: maxIndex, by: 1) {

                
                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append(((Double(i)*stepSize+xMin-(xMin+xMax)/2.0)*4.0*0.1))

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)

            }

            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
            
            
            
        case "Square Barrier":
            print("Getting Square Barrier")
        
            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
            let maxIndex = Int((xMax-stepSize)/stepSize)
            for i in stride(from: 1, through: Int(Double(maxIndex)*0.4), by: 1) {
        
                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append(0.0)
        
                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)
        
            }
        
            for i in stride(from: Int(Double(maxIndex)*0.4), to: Int(Double(maxIndex)*0.6), by: 1) {

                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append(15.000000001)

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)

            }

            for i in stride(from: Int(Double(maxIndex)*0.6), to: maxIndex, by: 1) {

                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append(0.0)

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)
            }
        
            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
            
            
            
        case "Triangle Barrier":
            print("Getting Triangle Barrier")
            var dataPoint: (xPoint: Double, yPoint: Double)
            var count = 0
            let maxIndex = Int((xMax-stepSize)/stepSize)
            
            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
        
            for i in stride(from: 1, to: Int(Double(maxIndex)*0.4), by: 1) {
//                print(i)
                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append(0.0)

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)
            }

            for i in stride(from: Int(Double(maxIndex)*0.4), to: Int(Double(maxIndex)*0.5), by: 1) {
//                print(i)
                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append((abs((Double(i)*stepSize+xMin)-(xMin + (xMax-xMin)*0.4))*3.0))

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)

            }

            for i in stride(from: Int(Double(maxIndex)*0.5), to: Int(Double(maxIndex)*0.6), by: 1) {
                print(i)
                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append((abs((Double(i)*stepSize+xMin)-(xMax - (xMax-xMin)*0.4))*3.0))

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)

            }

            for i in stride(from: Int(Double(maxIndex)*0.6), through: maxIndex, by: 1) {
//                print(i)
                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append(0.0)

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)
            }
            
            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
            
            
        case "Coupled Parabolic Well":
            print("Getting Coupled Parabolic Well")
            var dataPoint: (xPoint: Double, yPoint: Double)
            var count = 0

            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
            let maxIndex = Int((xMax-stepSize)/stepSize)

            for i in stride(from: 1, to: Int(Double(maxIndex)*0.5), by: 1) {

                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append((pow((Double(i)*stepSize+xMin-(xMin+(xMax-xMin)/4.0)), 2.0)))
//                print((pow((Double(i)*stepSize+xMin-(xMin+(xMax-xMin)/4.0)), 2.0)))

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)

            }

            for i in stride(from: Int(Double(maxIndex)*0.5), through: maxIndex, by: 1) {

                potential.xArray.append(Double(i)*stepSize+xMin)
                potential.yArray.append((pow(((Double(i)*stepSize+xMin)-(xMax-(xMax-xMin)/4.0)), 2.0)))
//                print((pow(((Double(i)*stepSize+xMin)-(xMax-(xMax-xMin)/4.0)), 2.0)))

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)

            }

            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
    
        case "Coupled Square Well + Field":
            print("Getting Coupled Square Well + Field")
            var dataPoint: (xPoint: Double, yPoint: Double)

            
            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
            let maxIndex = Int((xMax-stepSize)/stepSize)
            for i in stride(from: 1, to: Int(Double(maxIndex)*0.4), by: 1) {

                potential.xArray.append(Double(i)*stepSize + xMin)
                potential.yArray.append(0.0)

            }

            for i in stride(from: Int(Double(maxIndex)*0.4), to: Int(Double(maxIndex)*0.6), by: 1) {

                potential.xArray.append(Double(i)*stepSize + xMin)
                potential.yArray.append(4.0)

            }

            for i in stride(from: Int(Double(maxIndex)*0.6), to: maxIndex, by: 1) {

                potential.xArray.append(Double(i)*stepSize + xMin)
                potential.yArray.append(0.0)

            }

            for i in 1 ..< (potential.xArray.count) {

                potential.yArray[i] += ((potential.xArray[i]-xMin)*4.0*0.1)
                dataPoint = (xPoint: potential.xArray[i], yPoint: potential.yArray[i])
                contentArray.append(dataPoint)
            }

            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
                                        
                                        
        case "Harmonic Oscillator":
            print("Getting Harmonic Oscillator")
            var dataPoint: (xPoint: Double, yPoint: Double)
            var count = 0

            let xMinHO = -20.0
            let xMaxHO = 20.0
            let xStepHO = 0.001

            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)

            for i in stride(from: xMinHO+xStepHO, through: xMaxHO-xStepHO, by: xStepHO) {

                potential.xArray.append(i+xMaxHO)
                potential.yArray.append((pow((i-(xMaxHO+xMinHO)/2.0), 2.0)/15.0))

                count = potential.xArray.count
                dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                contentArray.append(dataPoint)
            }

            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
                                        
                       
        case "Kronig - Penney":
            print("Getting Kronig Penny")

            var dataPoint: (xPoint: Double, yPoint: Double)
            var count = 0
            let xMinKP = 0.0
            let xStepKP = 0.001
            let numberOfBarriers = 10.0
            let boxLength = 10.0
            let barrierPotential = 100.0*hbar2overm/2.0
            let latticeSpacing = boxLength/numberOfBarriers
            let barrierWidth = 1.0/6.0*latticeSpacing
            var barrierNumber = 1;
            var currentBarrierPosition = 0.0
            var inBarrier = false;
            let xMaxKP = boxLength


            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)

            for i in stride(from: xMinKP+xStepKP, through: xMaxKP-xStepKP, by: xStepKP) {

                currentBarrierPosition = -latticeSpacing/2.0 + Double(barrierNumber)*latticeSpacing

                if( (abs(i-currentBarrierPosition)) < (barrierWidth/2.0)) {

                    inBarrier = true

//                    potential.onedArray.append((xCoord: i, Potential: barrierPotential))

                    potential.xArray.append(i)
                    potential.yArray.append(barrierPotential)

                    count = potential.xArray.count
                    dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                    contentArray.append(dataPoint)


                }
                else {

                    if (inBarrier){

                        inBarrier = false
                        barrierNumber += 1

                    }

                    potential.xArray.append(i)
                    potential.yArray.append(0.0)

                    count = potential.xArray.count
                    dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
                    contentArray.append(dataPoint)
                }
            }

            potential.xArray.append(xMax)
            potential.yArray.append(5000000.0)

            dataPoint = (xPoint: potential.xArray[count-1], yPoint: potential.yArray[count-1])
            contentArray.append(dataPoint)
            
        default:
            //We will have the potential default be zero
            print("Getting Default")
            for i in stride(from: xMin, through: xMax, by: stepSize){
                potential.onedArray.removeAll()
                potential.xArray.append(i)
                potential.yArray.append(0.0)
            }
        }
    }
}




//func getPotential(potentialType: String, xMin: Double, xMax: Double, xStep: Double)
//    {
//        potential.oneDPotentialArray.removeAll()
//        potential.xArray.removeAll()
//        potential.yArray.removeAll()
//
//        xOffset = 0.0
//
//        var dataPoint: plotDataType = [:]
//        var count = 0
//
//        switch potentialType {
//        case "Square Well":
//
//                startPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//            for i in stride(from: xMin+xStep, through: xMax-xStep, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//            }
//
//                finishPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//
//        case "Linear Well":
//
//            startPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//            for i in stride(from: xMin+xStep, through: xMax-xStep, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append((i-xMin)*4.0*1.3)
//                //potential.yArray.append((i-xMin)*0.25)
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            finishPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//        case "Parabolic Well":
//
//            startPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//            for i in stride(from: xMin+xStep, through: xMax-xStep, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append((pow((i-(xMax+xMin)/2.0), 2.0)/1.0))
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            finishPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//        case "Square + Linear Well":
//
//            startPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//            for i in stride(from: xMin+xStep, to: (xMax+xMin)/2.0, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: (xMin+xMax)/2.0, through: xMax-xStep, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(((i-(xMin+xMax)/2.0)*4.0*0.1))
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            finishPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//
//        case "Square Barrier":
//
//            startPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//            for i in stride(from: xMin+xStep, to: xMin + (xMax-xMin)*0.4, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.4, to: xMin + (xMax-xMin)*0.6, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(15.000000001)
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.6, to: xMax, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//            }
//
//            finishPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//        case "Triangle Barrier":
//
//            var dataPoint: plotDataType = [:]
//            var count = 0
//
//            startPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//            for i in stride(from: xMin+xStep, to: xMin + (xMax-xMin)*0.4, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.4, to: xMin + (xMax-xMin)*0.5, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append((abs(i-(xMin + (xMax-xMin)*0.4))*3.0))
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.5, to: xMin + (xMax-xMin)*0.6, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append((abs(i-(xMax - (xMax-xMin)*0.4))*3.0))
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.6, to: xMax, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//            }
//
//            finishPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//        case "Coupled Parabolic Well":
//
//            var dataPoint: plotDataType = [:]
//            var count = 0
//
//            startPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//            for i in stride(from: xMin+xStep, to: xMin + (xMax-xMin)*0.5, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append((pow((i-(xMin+(xMax-xMin)/4.0)), 2.0)))
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.5, through: xMax-xStep, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append((pow((i-(xMax-(xMax-xMin)/4.0)), 2.0)))
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            finishPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//        case "Coupled Square Well + Field":
//
//            var dataPoint: plotDataType = [:]
//
//            startPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//            for i in stride(from: xMin+xStep, to: xMin + (xMax-xMin)*0.4, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.4, to: xMin + (xMax-xMin)*0.6, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(4.0)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.6, to: xMax, by: xStep) {
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//            }
//
//            for i in 1 ..< (potential.xArray.count) {
//
//                potential.yArray[i] += ((potential.xArray[i]-xMin)*4.0*0.1)
//                dataPoint = [.X: potential.xArray[i], .Y: potential.yArray[i]]
//                contentArray.append(dataPoint)
//            }
//
//
//            finishPotential(xMin: xMin, xMax: xMax, xStep: xStep)
//
//        case "Harmonic Oscillator":
//
//            var dataPoint: plotDataType = [:]
//            var count = 0
//
//            let xMinHO = -20.0
//            let xMaxHO = 20.0
//            let xStepHO = 0.001
//
//            startPotential(xMin: xMinHO+xMaxHO, xMax: xMaxHO+xMaxHO, xStep: xStepHO)
//
//            for i in stride(from: xMinHO+xStepHO, through: xMaxHO-xStepHO, by: xStepHO) {
//
//                potential.xArray.append(i+xMaxHO)
//                potential.yArray.append((pow((i-(xMaxHO+xMinHO)/2.0), 2.0)/15.0))
//
//                count = potential.xArray.count
//                dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                contentArray.append(dataPoint)
//            }
//
//            finishPotential(xMin: xMinHO+xMaxHO, xMax: xMaxHO+xMaxHO, xStep: xStepHO)
//
//        case "Kronig - Penney":
//
//            var dataPoint: plotDataType = [:]
//            var count = 0
//
//            let xMinKP = 0.0
//
//            let xStepKP = 0.001
//
//            let numberOfBarriers = 10.0
//            let boxLength = 10.0
//            let barrierPotential = 100.0*hbar2overm/2.0
//            let latticeSpacing = boxLength/numberOfBarriers
//            let barrierWidth = 1.0/6.0*latticeSpacing
//            var barrierNumber = 1;
//            var currentBarrierPosition = 0.0
//            var inBarrier = false;
//
//            let xMaxKP = boxLength
//
//
//            startPotential(xMin: xMinKP, xMax: xMaxKP, xStep: xStepKP)
//
//            for i in stride(from: xMinKP+xStepKP, through: xMaxKP-xStepKP, by: xStepKP) {
//
//                currentBarrierPosition = -latticeSpacing/2.0 + Double(barrierNumber)*latticeSpacing
//
//                if( (abs(i-currentBarrierPosition)) < (barrierWidth/2.0)) {
//
//                    inBarrier = true
//
//                    potential.oneDPotentialArray.append((xCoord: i, Potential: barrierPotential))
//
//                    potential.xArray.append(i)
//                    potential.yArray.append(barrierPotential)
//
//                    count = potential.xArray.count
//                    dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                    contentArray.append(dataPoint)
//
//
//                }
//                else {
//
//                    if (inBarrier){
//
//                        inBarrier = false
//                        barrierNumber += 1
//
//                    }
//
//                    potential.xArray.append(i)
//                    potential.yArray.append(0.0)
//
//                    count = potential.xArray.count
//                    dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//                    contentArray.append(dataPoint)
//
//
//                }
//
//
//            }
//
//            potential.xArray.append(xMax)
//            potential.yArray.append(5000000.0)
//
//            dataPoint = [.X: potential.xArray[count-1], .Y: potential.yArray[count-1]]
//            contentArray.append(dataPoint)
//
//            /** Fixes Bug In Plotting Library not displaying the last point **/
//            dataPoint = [.X: xMax+xStep, .Y: 5000000.0]
//            contentArray.append(dataPoint)
//
//            let xMin = potential.minX(minArray: potential.xArray)
//            let xMax = potential.maxX(maxArray: potential.xArray)
//            let yMin = potential.minY(minArray: potential.yArray)
//            var yMax = potential.maxY(maxArray: potential.yArray)
//
//            if yMax > 500 { yMax = 10}
//
//            makePlot(xLabel: "x Å", yLabel: "Potential V", xMin: (xMin - 1.0), xMax: (xMax + 1.0), yMin: yMin-1.2, yMax: yMax+0.2)
//
//            contentArray.removeAll()
//
//        case "Variable Kronig - Penney":
//
//            /****  Get Parameters ****/
//
//            if let kpDataController = storyboard!.instantiateController(withIdentifier: "theSecondViewController") as? secondViewController {
//                kpDataController.delegate = self
//                presentAsSheet(kpDataController)
//            }
//
//
//        case "KP2-a":
//
//            var dataPoint: plotDataType = [:]
//            var count = 0
//
//            let xMinKP = 0.0
//
//            let xStepKP = 0.001
//
//           // let numberOfBarriers = 2.0
//            let boxLength = 10.0
//            let barrierPotential = 100.0*hbar2overm/2.0
//            let latticeSpacing = 1.0 //boxLength/numberOfBarriers
//            let barrierWidth = 1.0/6.0*latticeSpacing
//            var barrierNumber = 1;
//            var currentBarrierPosition = 0.0
//            var inBarrier = false;
//
//            let xManKP = boxLength
//
//
//            potential.oneDPotentialArray.append((xCoord: xMinKP, Potential: 5000000.0))
//            dataPoint = [.X: potential.oneDPotentialArray[0].xCoord, .Y: potential.oneDPotentialArray[0].Potential]
//            contentArray.append(dataPoint)
//
//            for i in stride(from: xMinKP+xStepKP, through: xManKP-xStepKP, by: xStepKP) {
//
//                let term = (-latticeSpacing/2.0) * (pow(-1.0, Double(barrierNumber))) - Double(barrierNumber)*Double(barrierNumber-1) * (pow(-1.0, Double(barrierNumber)))
//
//                currentBarrierPosition =  term + Double(barrierNumber)*latticeSpacing*4.0
//
//                if( (abs(i-currentBarrierPosition)) < (barrierWidth/2.0)) {
//
//                    inBarrier = true
//
//                    potential.oneDPotentialArray.append((xCoord: i, Potential: barrierPotential))
//
//                    let count = potential.oneDPotentialArray.count - 1
//                    let dataPoint: plotDataType = [.X: potential.oneDPotentialArray[count].xCoord, .Y: potential.oneDPotentialArray[count].Potential]
//                    contentArray.append(dataPoint)
//
//                }
//                else {
//
//                    if (inBarrier){
//
//                        inBarrier = false
//                        barrierNumber += 1
//
//                    }
//
//                    potential.oneDPotentialArray.append((xCoord: i, Potential: 0.0))
//
//                    let count = potential.oneDPotentialArray.count - 1
//                    let dataPoint: plotDataType = [.X: potential.oneDPotentialArray[count].xCoord, .Y: potential.oneDPotentialArray[count].Potential]
//                    contentArray.append(dataPoint)
//
//
//                }
//
//
//            }
//
//            count = potential.oneDPotentialArray.count
//            potential.oneDPotentialArray.append((xCoord: xManKP, Potential: 5000000.0))
//            dataPoint = [.X: potential.oneDPotentialArray[count-1].xCoord, .Y: potential.oneDPotentialArray[count-1].Potential]
//            contentArray.append(dataPoint)
//
//            /** Fixes Bug In Plotting Library not displaying the last point **/
//            dataPoint = [.X: xManKP+xStepKP, .Y: 5000000]
//            contentArray.append(dataPoint)
//
//            let xMin = potential.minX(minArray: potential.xArray)
//            let xMax = potential.maxX(maxArray: potential.xArray)
//            let yMin = potential.minY(minArray: potential.yArray)
//            var yMax = potential.maxY(maxArray: potential.yArray)
//
//            if yMax > 500 { yMax = 10}
//
//            makePlot(xLabel: "x Å", yLabel: "Potential V", xMin: (xMin - 1.0), xMax: (xMax + 1.0), yMin: yMin-1.2, yMax: yMax+0.2)
//
//            contentArray.removeAll()
//
//        default:
//            let tab: Character = "\t"
//            let geFilePanel: NSOpenPanel = NSOpenPanel()
//            var filePath :URL = URL(string:("file://"))!
//
//            var dataPoint: plotDataType = [:]
//
//            geFilePanel.runModal()
//
//            // Get the file path from the NSSavePanel
//
//            filePath = URL(string:("file://" + (geFilePanel.url?.path)!))!
//
//            print(filePath)
//
//            do {
//                let tsv = try CSV(url: filePath, delimiter: tab, encoding: String.Encoding.utf8, loadColumns: true)
//
//                var xArray: Array = (tsv.namedColumns[tsv.header[0]] as Array?)!
//                var yArray: Array = (tsv.namedColumns[tsv.header[1]] as Array?)!
//
//
//                for index in 0..<xArray.count {
//
//                    potential.xArray.append(Double(xArray[index])!)
//                    potential.yArray.append(Double(yArray[index])!)
//
//                }
//
//                let xMin = potential.minX(minArray: potential.xArray)
//                let xMax = potential.maxX(maxArray: potential.xArray)
//                let yMin = potential.minY(minArray: potential.yArray)
//                var yMax = potential.maxY(maxArray: potential.yArray)
//
//                if (xMin < 0.0) {
//
//                    xOffset = -xMin
//
//                    for i in 0..<potential.xArray.count {
//
//                        dataPoint = [.X: potential.xArray[i], .Y: potential.yArray[i]]
//                        contentArray.append(dataPoint)
//
//                        potential.xArray[i] += xOffset
//
//                    }
//
//
//                }
//
//                if yMax > 500 { yMax = 10}
//
//                makePlot(xLabel: "x Å", yLabel: "Potential V", xMin: (xMin - 1.0), xMax: (xMax + 1.0), yMin: yMin-1.2, yMax: yMax+0.2)
//
//                contentArray.removeAll()
//
//
//            } catch {
//                // Error handling
//            }
//
//
//
//        }
//
//
//    }
//
//    func userDidEnterInformation(info: String?) {
//
//        setupVariableKronigPenney(numberOfBarriers: Double(info!)!)
//
//    }
//
//
//    func setupVariableKronigPenney(numberOfBarriers: Double) {
//        var dataPoint: plotDataType = [:]
//        var count = 0
//
//        let xMinKP = 0.0
//
//        let xStepKP = 0.001
//
//        let boxLength = numberOfBarriers
//        let barrierPotential = 100.0*hbar2overm/2.0
//        let latticeSpacing = boxLength/numberOfBarriers
//        let barrierWidth = 1.0/6.0*latticeSpacing
//        var barrierNumber = 1;
//        var currentBarrierPosition = 0.0
//        var inBarrier = false;
//
//        let xMaxKP = boxLength
//
//
//        startPotential(xMin: xMinKP, xMax: xMaxKP, xStep: xStepKP)
//
//        for i in stride(from: xMinKP+xStepKP, through: xMaxKP-xStepKP, by: xStepKP) {
//
//            currentBarrierPosition = -latticeSpacing/2.0 + Double(barrierNumber)*latticeSpacing
//
//            if( (abs(i-currentBarrierPosition)) < (barrierWidth/2.0)) {
//
//                inBarrier = true
//
//                potential.xArray.append(i)
//                potential.yArray.append(barrierPotential)
//
//                count = potential.xArray.count - 1
//                dataPoint = [.X: potential.xArray[count], .Y: potential.yArray[count]]
//                contentArray.append(dataPoint)
//
//            }
//            else {
//
//                if (inBarrier){
//
//                    inBarrier = false
//                    barrierNumber += 1
//
//                }
//
//                potential.xArray.append(i)
//                potential.yArray.append(0.0)
//
//                count = potential.xArray.count - 1
//                dataPoint = [.X: potential.xArray[count], .Y: potential.yArray[count]]
//                contentArray.append(dataPoint)
//
//
//            }
//
//
//        }
//
//        count = potential.xArray.count-1
//        potential.xArray.append(xMaxKP)
//        potential.yArray.append(5000000.0)
//        dataPoint = [.X: potential.xArray[count], .Y: potential.yArray[count]]
//        contentArray.append(dataPoint)
//
//        /** Fixes Bug In Plotting Library not displaying the last point **/
//        dataPoint = [.X: xMaxKP+xStepKP, .Y: 5000000]
//        contentArray.append(dataPoint)
//
//        let xMin = potential.minX(minArray: potential.xArray)
//        let xMax = potential.maxX(maxArray: potential.xArray)
//        let yMin = potential.minY(minArray: potential.yArray)
//        var yMax = potential.maxY(maxArray: potential.yArray)
//
//        if yMax > 500 { yMax = 10}
//
//        makePlot(xLabel: "x Å", yLabel: "Potential V", xMin: (xMin - 1.0), xMax: (xMax + 1.0), yMin: yMin-1.2, yMax: yMax+0.2)
//
//        contentArray.removeAll()
//
//    }
