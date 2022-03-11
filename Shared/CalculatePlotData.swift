//
//  CalculatePlotData.swift
//  SwiftUICorePlotExample
//
//  Created by Jeff Terry on 12/22/20.
//

import Foundation
import SwiftUI
import CorePlot

class CalculatePlotData: ObservableObject {
    
    var plotDataModel: PlotDataClass? = nil
    var theText = ""
    

    @MainActor func setThePlotParameters(color: String, xLabel: String, yLabel: String, title: String, minX: Double, maxX: Double, maxY: Double) {
        //set the Plot Parameters
        plotDataModel!.changingPlotParameters.yMax = maxY+1.0
        plotDataModel!.changingPlotParameters.yMin = -1.0 * maxY - 1.0
        plotDataModel!.changingPlotParameters.xMax = maxX
        plotDataModel!.changingPlotParameters.xMin = minX
        plotDataModel!.changingPlotParameters.xLabel = xLabel
        plotDataModel!.changingPlotParameters.yLabel = yLabel
        
        if color == "Red"{
            plotDataModel!.changingPlotParameters.lineColor = .red()
        }
        else{
            
            plotDataModel!.changingPlotParameters.lineColor = .blue()
        }
        plotDataModel!.changingPlotParameters.title = title
        
        plotDataModel!.zeroData()
    }
    
    @MainActor func appendDataToPlot(plotData: [plotDataType]) {
        plotDataModel!.appendData(dataPoint: plotData)
    }
    
    func plotYEqualsX() async
    {
        
        theText = "y = x\n"
        
        await setThePlotParameters(color: "Red", xLabel: "x", yLabel: "y", title: "y = x", minX: 0.0, maxX: 10.0, maxY: 10.0)
        
        await resetCalculatedTextOnMainThread()
        
        
        var plotData :[plotDataType] =  []
        
        
        for i in 0 ..< 120 {
             
            //create x values here

            let x = -2.0 + Double(i) * 0.2

        //create y values here

        let y = x


            let dataPoint: plotDataType = [.X: x, .Y: y]
            plotData.append(contentsOf: [dataPoint])
            theText += "x = \(x), y = \(y)\n"
        
        }
        
        await appendDataToPlot(plotData: plotData)
        await updateCalculatedTextOnMainThread(theText: theText)
        
        
    }
    
    
    func plotFunction(data: [(xPoint: Double, yPoint: Double)], labelForX: String, titleStr: String, lowerE: Double, upperE: Double) async
    {
        
        //set the Plot Parameters
        var ymax = 0.0
        
        for i in data{
            if i.yPoint > ymax{
                ymax = i.yPoint
            }
            
        }
        print(ymax)
        
        await plotDataModel!.changingPlotParameters.yMax = ymax+1.0
        await plotDataModel!.changingPlotParameters.yMin = -1.0 * ymax - 1.0
        await plotDataModel!.changingPlotParameters.xMax = upperE+1.0
        await plotDataModel!.changingPlotParameters.xMin = lowerE-1.0
        await plotDataModel!.changingPlotParameters.xLabel = labelForX
        await plotDataModel!.changingPlotParameters.yLabel = "Psi"
        await plotDataModel!.changingPlotParameters.lineColor = .blue()
        await plotDataModel!.changingPlotParameters.title = titleStr
        await plotDataModel!.zeroData()
        
        await setThePlotParameters(color: "Blue", xLabel: labelForX, yLabel: "Psi", title: titleStr, minX: lowerE-1.0, maxX: upperE+1.0, maxY: ymax)
        
        await resetCalculatedTextOnMainThread()
        
        theText = ""
        
        var plotData :[plotDataType] =  []
        for i in data{
            
            
        //create x values here
        let x = i.xPoint

        //create y values here
        let y = i.yPoint
        
        let dataPoint: plotDataType = [.X: x, .Y: y]
        plotData.append(contentsOf: [dataPoint])
        theText += "x = \(x), y = \(y)\n"
            
        }
        
        await appendDataToPlot(plotData: plotData)
        await updateCalculatedTextOnMainThread(theText: theText)
        
        return
    }
    
    
        @MainActor func resetCalculatedTextOnMainThread() {
            //Print Header
            plotDataModel!.calculatedText = ""
    
        }
    
    
        @MainActor func updateCalculatedTextOnMainThread(theText: String) {
            //Print Header
            plotDataModel!.calculatedText += theText
        }
    
}



