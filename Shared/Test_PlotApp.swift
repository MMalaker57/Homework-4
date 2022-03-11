//
//  Test_PlotApp.swift
//  Shared
//
//  Created by Jeff Terry on 1/25/21.
//

import SwiftUI

@main
struct Test_PlotApp: App {
    
    @StateObject var plotData = PlotClass()
    
    var body: some Scene {
        WindowGroup {
            TabView {
                ContentView(functional: [(xPoint: 0.0, yPoint: 0.0)])
                    .environmentObject(plotData)
                    .tabItem {
                        Text("Plot")
                    }
//                TextView()
//                    .environmentObject(plotData)
//                    .tabItem {
//                        Text("Text")
//                    }
//
                            
            }
            
        }
    }
}
