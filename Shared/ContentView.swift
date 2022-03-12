//
//  ContentView.swift
//  Shared
//
//  Created by Jeff Terry on 1/25/21.
//

import SwiftUI
import CorePlot

typealias plotDataType = [CPTScatterPlotField : Double]

struct ContentView: View {
    @EnvironmentObject var plotData :PlotClass
    @ObservedObject var plotDataModel = PlotDataClass(fromLine: true)
    
    @ObservedObject private var calculator = CalculatePlotData()
    @ObservedObject var schrodinger = Schrodinger()
    @ObservedObject var potentialV = Potential()
    @State var isChecked:Bool = false
    @State var tempInput = ""
    @State var selector = 0
    
    @State var lowerX: Double = 0.0
    @State var upperX = 10.0
    @State var maxE = 30.0
    @State var eStep = 0.1
    @State var stepInX = 0.1
    
    @State var functional: [(xPoint: Double, yPoint: Double)]
    @State var functionalData: (E: [Double],psi: [Double]) = (E: [0.0],psi: [0.0])
    
    @State var lowerXStr = ""
    @State var upperXStr = ""
    @State var maxEStr = ""
    @State var eStepStr = ""
    @State var xStepStr = ""
    @State var VString = ""
    @State var eigenString = ""
    @State var VToUse = ""
    @State var potentials: [String] = ["Square Well", "Linear Well","Parabolic Well","Square + Linear Well","Square Barrier","Triangle Barrier","Coupled Parabolic Well","Harmonic Oscillator","Kronig - Penney"]
            
    
    @State var energies: [String] = ["0 eV"]
    
    @State var psiArray: [(E: Double, psi: [Double])] = []
    @State var selectedPsi: [(E: Double, psi: [Double])] = []

    var body: some View {
        
        VStack{
      
            CorePlot(dataForPlot: $plotData.plotArray[selector].plotData, changingPlotParameters: $plotData.plotArray[selector].changingPlotParameters)
                .setPlotPadding(left: 10)
                .setPlotPadding(right: 10)
                .setPlotPadding(top: 10)
                .setPlotPadding(bottom: 10)
                .padding()
            
            Divider()
            /*
                     CorePlot(dataForPlot: $plotData.plotArray[0].plotData, changingPlotParameters: $plotData.plotArray[0].changingPlotParameters)
                         .setPlotPadding(left: 10)
                         .setPlotPadding(right: 10)
                         .setPlotPadding(top: 10)
                         .setPlotPadding(bottom: 10)
                         .padding()
                     
                     Divider()
                     */
   /*
            CorePlot(dataForPlot: $plotData.plotArray[1].plotData, changingPlotParameters: $plotData.plotArray[1].changingPlotParameters)
                .setPlotPadding(left: 10)
                .setPlotPadding(right: 10)
                .setPlotPadding(top: 10)
                .setPlotPadding(bottom: 10)
                .padding()
            
            Divider()
            */
            VStack{
                HStack{
                    Text("Input Min X")
                    TextField("", text: $lowerXStr)
                        .padding()
                }
                HStack{
                    Text("Input Max X")
                    TextField("", text: $upperXStr)
                        .padding()
                }
                HStack{
                    Text("Input Max E")
                    TextField("", text: $maxEStr)
                        .padding()
                }
                HStack{
                    Text("Input E Step")
                    TextField("", text: $eStepStr)
                        .padding()
                }
                HStack{
                    Text("Input X Step")
                    TextField("", text: $xStepStr)
                        .padding()
                }
                HStack{
                    Picker("Potential", selection: $VString){
                        ForEach(potentials, id: \.self){
                            Text("\($0)").tag("\($0)")
                        }
                    }
                    
                    Picker("Energy", selection: $eigenString){
                        ForEach(energies, id: \.self){
                            Text("\($0)").tag("\($0)")
                        }
                    }
                    
                }
                .padding()
                
            }
            
            
            HStack{
                Button("Calculate", action: {
                    
                    self.setParams(lowerXStr: lowerXStr, upperXStr: upperXStr, eMaxStr: maxEStr, eStepStr: eStepStr, xStepStr: xStepStr)
                    energies.removeAll()
                    
                    self.setV(choice: VToUse)
                    schrodinger.setV(choice: eigenString, xMin: lowerX, xMax: upperX, stepSize: stepInX)
                    Task.init{ psiArray = await self.calculatepsiForAllE(minX: lowerX, maxX: upperX, maxE: maxE, xStep: stepInX, eStep: eStep)
                    }
                    
                    
                    
                }
                
                
                
                )
                .padding()
                
                Button("Draw", action: {
                    
                    Task.init{
                    await self.setPlotData(choice: eigenString)
                    }
                    self.plotData.objectWillChange.send()
                }
                )
                    .padding()
                Button("Clear", action: {
                    setupPlotDataModel(selector: 0)
                    Task.init{
                        await calculator.plotFunction(data: [(xPoint: 0.1, yPoint: 0.1)], labelForX: "x", titleStr: "Psi Vs X", lowerE: lowerX, upperE: upperX)
                    }
                    self.plotData.objectWillChange.send()
                }
                )
                    .padding()
                
            }
            
            
            
            
        }
        
        
        
    }
    
    @MainActor func setupPlotDataModel(selector: Int){
        
        calculator.plotDataModel = self.plotData.plotArray[selector]
    }
    
    
    
    
    func calculatepsiForAllE(minX: Double, maxX: Double, maxE: Double, xStep: Double, eStep: Double) async->[(E: Double, psi: [Double])]{
        functional.removeAll()
        potentialV.getPotential(potentialToGet: VToUse, xMin: minX, xMax: maxX, stepSize: xStep)
        schrodinger.setV(choice: VToUse, xMin: minX, xMax: maxX, stepSize: xStep)
        let psiStruct = await withTaskGroup(of: (E: Double, psi: [Double]).self, returning: [(E: Double, psi: [Double])].self, body: {taskGroup in
            
            var psiData: [(E: Double, psi: [Double])] = []
            for E in stride(from: 0+eStep, through: maxE, by: eStep){
                taskGroup.addTask{
//                    print(E)
                    let psiTemp = await schrodinger.calculatePsiForE(xMin: minX, xMax: maxX, stepSize: xStep, potentialStr: VString, E: E)
                    
                    return((E: E, psi: psiTemp))
                    
                }
            }
            
            var interimResults: [(E: Double, psi: [Double])] = []
            for await result in taskGroup{
                interimResults.append(result)
            }
            
            psiData = interimResults.sorted(by: {$0.E < $1.E})
//            for i in psiData{
//                print(i.E)
//                print(Double(i.psi.last ?? 50.0))
//            }
            return psiData
            
        })
        
        print("psiStruct is \(psiStruct.count) long")
        
        for i in psiStruct{
            functional.append((xPoint: i.E, yPoint: Double(i.psi.last ?? 50.0)))
        }
        
//        print(functional)
        
        let psiStruct2 = schrodinger.separateGoodE(minX: minX, maxX: maxX, stepSize: xStep, potentialStr: VString, data: psiStruct)
//        print(psiStruct2)
        energies.removeAll()
        energies.append("Functional")
        energies.append("Potential")
//        print(psiStruct2)
        
        
        var psiStruct3: [(E: Double, psi: [Double])] = []
        for i in stride(from: 0, to: psiStruct2.count, by: 1){
            energies.append("\(psiStruct2[i].E) eV")
            psiStruct3.append(schrodinger.normalize(data: psiStruct2[i], stepSize: xStep))
            
//            print(i.E)
        }
        
        return psiStruct3
    }
    
    func setParams(lowerXStr: String, upperXStr: String, eMaxStr: String, eStepStr: String, xStepStr: String){
     
        lowerX = Double(lowerXStr) ?? 0.0
        upperX = Double(upperXStr) ?? 1.0
        maxE = Double(eMaxStr) ?? 30.0
        eStep = Double(eStepStr) ?? 0.1
        stepInX = Double(xStepStr) ?? 0.1
        
        
    }
    
    func setPlotData(choice: String) async{
        setupPlotDataModel(selector: 0)
        
        print(eigenString)
        var plotData: [(xPoint: Double, yPoint: Double)] = []
        var count = 0
        var pickedPsi: [Double] = []
        await calculator.plotFunction(data: plotData, labelForX: "x", titleStr: "Psi vs X", lowerE: lowerX, upperE: upperX)
        
        
        if eigenString == "Functional"{
            
            plotData = functional
//            print(plotData)
            await calculator.plotFunction(data: plotData, labelForX: "E", titleStr: "Psi At A vs Energy", lowerE: 0.0, upperE: maxE)
            self.plotData.objectWillChange.send()
        }
        else{
            if eigenString == "Potential"{
                plotData = schrodinger.potential.contentArray
                await calculator.plotFunction(data: plotData, labelForX: "X", titleStr: "Potential vs X", lowerE: lowerX, upperE: upperX)
                self.plotData.objectWillChange.send()
            }
            else{
        for i in energies{
            
            if i == "0 eV"{
                pickedPsi = [0.0]
            }
            if i == eigenString{

                pickedPsi = psiArray[count-2].psi
//                print(pickedPsi)
                count = 0
                break
            }
            count += 1
            
        }
            
            count=0
        for y in pickedPsi{
//            let dataPoint: plotDataType = [.X: (lowerX + Double(count)*stepInX), .Y: y]
            plotData.append((xPoint: (lowerX + Double(count)*stepInX), yPoint: y))
//            data.append((xPoint: (lowerX + Double(count)*stepInX), yPoint: y))
            count+=1
        }
            count=0
            await calculator.plotFunction(data: plotData, labelForX: "x", titleStr: "Psi vs X", lowerE: lowerX, upperE: upperX)
            self.plotData.objectWillChange.send()
        
        }
        }
    }
    
    
    func setV(choice: String){
        print(VString)
        switch VString{
        case "Square Well":
            VToUse = "Square Well"
            print(VToUse)
        case "Linear Well":
            VToUse = "Linear Well"
            print(VToUse)
            
        case "Parabolic Well":
            VToUse = "Parabolic Well"
            print(VToUse)
            
        case "Square + Linear Well":
            VToUse = "Square + Linear Well"
            print(VToUse)
            
        case "Square Barrier":
            VToUse = "Square Barrier"
            print(VToUse)
            
        case "Triangle Barrier":
            VToUse = "Triangle Barrier"
            print(VToUse)
            
        case "Coupled Parabolic Well":
            VToUse = "Coupled Parabolic Well"
            print(VToUse)
            
        case "Harmonic Oscillator":
            VToUse = "Harmonic Oscillator"
            print(VToUse)
            
        case "Kronig - Penney":
            VToUse = "Kronig - Penney"
            print(VToUse)
        default:
            print("Resorting to Square Well")
            VToUse = "Square Well"
        }
    }
    
    /// calculate
    /// Function accepts the command to start the calculation from the GUI
//    func calculate() async {
//
//        //pass the plotDataModel to the Calculator
//       // calculator.plotDataModel = self.plotData.plotArray[0]
//
//        setupPlotDataModel(selector: 0)
//
//     //   Task{
//
//
//            let _ = await withTaskGroup(of:  Void.self) { taskGroup in
//
//
//
//                taskGroup.addTask {
//
//
//        var temp = 0.0
//
//
//
//        //Calculate the new plotting data and place in the plotDataModel
//        await calculator.plotFunction(data: <#T##[(xPoint: Double, yPoint: Double)]#>)()
//
//                    // This forces a SwiftUI update. Force a SwiftUI update.
//        await self.plotData.objectWillChange.send()
//
//                }
//
//
//            }
//
//  //      }
//
//
//    }
    
    /// calculate
    /// Function accepts the command to start the calculation from the GUI
//    func calculate2() async {
//
//
//        //pass the plotDataModel to the Calculator
//       // calculator.plotDataModel = self.plotData.plotArray[0]
//
//        setupPlotDataModel(selector: 1)
//
//     //   Task{
//
//
//            let _ = await withTaskGroup(of:  Void.self) { taskGroup in
//
//
//
//                taskGroup.addTask {
//
//
//        var temp = 0.0
//
//
//
//        //Calculate the new plotting data and place in the plotDataModel
//        await calculator.plotYEqualsX()
//
//                    // This forces a SwiftUI update. Force a SwiftUI update.
//        await self.plotData.objectWillChange.send()
//
//                }
//
//            }
//
//    //    }
//
//
//
//    }
    

   
    
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView( functional: [(xPoint: 0.0, yPoint: 0.0)])
    }
}
