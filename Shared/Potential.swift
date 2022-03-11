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
    
    
    func getPotential(potentialToGet: String, xMin: Double, xMax: Double, stepSize: Double){
        
        potential.onedArray = []
        potential.xArray = []
        potential.yArray = []
        
        
        
        switch potentialToGet{
            
        case "Square Well":
            //We set the bounds of the potential as "infinite" and the rest as zero
            potential.xArray.append(xMin)
            potential.yArray.append(1000000000000000)
            //Potential is in index space, not x space
            let maxIndex = Int((xMax-stepSize)/stepSize)
            for i in stride(from: 0, through: maxIndex-1, by: 1){
                print(i)
                potential.xArray.append(Double(i))
                potential.yArray.append(0.0)
            }
            potential.xArray.append(xMax)
            potential.yArray.append(1000000000000000)
            
            
        default:
            //We will have the potential default be zero
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
//        potential.oneDPotentialXArray.removeAll()
//        potential.oneDPotentialYArray.removeAll()
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append((i-xMin)*4.0*1.3)
//                //potential.oneDPotentialYArray.append((i-xMin)*0.25)
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append((pow((i-(xMax+xMin)/2.0), 2.0)/1.0))
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: (xMin+xMax)/2.0, through: xMax-xStep, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(((i-(xMin+xMax)/2.0)*4.0*0.1))
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.4, to: xMin + (xMax-xMin)*0.6, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(15.000000001)
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.6, to: xMax, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//                contentArray.append(dataPoint)
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.4, to: xMin + (xMax-xMin)*0.5, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append((abs(i-(xMin + (xMax-xMin)*0.4))*3.0))
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.5, to: xMin + (xMax-xMin)*0.6, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append((abs(i-(xMax - (xMax-xMin)*0.4))*3.0))
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.6, to: xMax, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append((pow((i-(xMin+(xMax-xMin)/4.0)), 2.0)))
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//                contentArray.append(dataPoint)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.5, through: xMax-xStep, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append((pow((i-(xMax-(xMax-xMin)/4.0)), 2.0)))
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.4, to: xMin + (xMax-xMin)*0.6, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(4.0)
//
//            }
//
//            for i in stride(from: xMin + (xMax-xMin)*0.6, to: xMax, by: xStep) {
//
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//            }
//
//            for i in 1 ..< (potential.oneDPotentialXArray.count) {
//
//                potential.oneDPotentialYArray[i] += ((potential.oneDPotentialXArray[i]-xMin)*4.0*0.1)
//                dataPoint = [.X: potential.oneDPotentialXArray[i], .Y: potential.oneDPotentialYArray[i]]
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
//                potential.oneDPotentialXArray.append(i+xMaxHO)
//                potential.oneDPotentialYArray.append((pow((i-(xMaxHO+xMinHO)/2.0), 2.0)/15.0))
//
//                count = potential.oneDPotentialXArray.count
//                dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                    potential.oneDPotentialXArray.append(i)
//                    potential.oneDPotentialYArray.append(barrierPotential)
//
//                    count = potential.oneDPotentialXArray.count
//                    dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
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
//                    potential.oneDPotentialXArray.append(i)
//                    potential.oneDPotentialYArray.append(0.0)
//
//                    count = potential.oneDPotentialXArray.count
//                    dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//                    contentArray.append(dataPoint)
//
//
//                }
//
//
//            }
//
//            potential.oneDPotentialXArray.append(xMax)
//            potential.oneDPotentialYArray.append(5000000.0)
//
//            dataPoint = [.X: potential.oneDPotentialXArray[count-1], .Y: potential.oneDPotentialYArray[count-1]]
//            contentArray.append(dataPoint)
//
//            /** Fixes Bug In Plotting Library not displaying the last point **/
//            dataPoint = [.X: xMax+xStep, .Y: 5000000.0]
//            contentArray.append(dataPoint)
//
//            let xMin = potential.minX(minArray: potential.oneDPotentialXArray)
//            let xMax = potential.maxX(maxArray: potential.oneDPotentialXArray)
//            let yMin = potential.minY(minArray: potential.oneDPotentialYArray)
//            var yMax = potential.maxY(maxArray: potential.oneDPotentialYArray)
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
//            let xMin = potential.minX(minArray: potential.oneDPotentialXArray)
//            let xMax = potential.maxX(maxArray: potential.oneDPotentialXArray)
//            let yMin = potential.minY(minArray: potential.oneDPotentialYArray)
//            var yMax = potential.maxY(maxArray: potential.oneDPotentialYArray)
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
//                    potential.oneDPotentialXArray.append(Double(xArray[index])!)
//                    potential.oneDPotentialYArray.append(Double(yArray[index])!)
//
//                }
//
//                let xMin = potential.minX(minArray: potential.oneDPotentialXArray)
//                let xMax = potential.maxX(maxArray: potential.oneDPotentialXArray)
//                let yMin = potential.minY(minArray: potential.oneDPotentialYArray)
//                var yMax = potential.maxY(maxArray: potential.oneDPotentialYArray)
//
//                if (xMin < 0.0) {
//
//                    xOffset = -xMin
//
//                    for i in 0..<potential.oneDPotentialXArray.count {
//
//                        dataPoint = [.X: potential.oneDPotentialXArray[i], .Y: potential.oneDPotentialYArray[i]]
//                        contentArray.append(dataPoint)
//
//                        potential.oneDPotentialXArray[i] += xOffset
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(barrierPotential)
//
//                count = potential.oneDPotentialXArray.count - 1
//                dataPoint = [.X: potential.oneDPotentialXArray[count], .Y: potential.oneDPotentialYArray[count]]
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
//                potential.oneDPotentialXArray.append(i)
//                potential.oneDPotentialYArray.append(0.0)
//
//                count = potential.oneDPotentialXArray.count - 1
//                dataPoint = [.X: potential.oneDPotentialXArray[count], .Y: potential.oneDPotentialYArray[count]]
//                contentArray.append(dataPoint)
//
//
//            }
//
//
//        }
//
//        count = potential.oneDPotentialXArray.count-1
//        potential.oneDPotentialXArray.append(xMaxKP)
//        potential.oneDPotentialYArray.append(5000000.0)
//        dataPoint = [.X: potential.oneDPotentialXArray[count], .Y: potential.oneDPotentialYArray[count]]
//        contentArray.append(dataPoint)
//
//        /** Fixes Bug In Plotting Library not displaying the last point **/
//        dataPoint = [.X: xMaxKP+xStepKP, .Y: 5000000]
//        contentArray.append(dataPoint)
//
//        let xMin = potential.minX(minArray: potential.oneDPotentialXArray)
//        let xMax = potential.maxX(maxArray: potential.oneDPotentialXArray)
//        let yMin = potential.minY(minArray: potential.oneDPotentialYArray)
//        var yMax = potential.maxY(maxArray: potential.oneDPotentialYArray)
//
//        if yMax > 500 { yMax = 10}
//
//        makePlot(xLabel: "x Å", yLabel: "Potential V", xMin: (xMin - 1.0), xMax: (xMax + 1.0), yMin: yMin-1.2, yMax: yMax+0.2)
//
//        contentArray.removeAll()
//
//    }
