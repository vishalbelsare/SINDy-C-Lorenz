function Error = SimulationError (Y_Actual, Y_Simulated)

    Error = (Y_Simulated - Y_Actual) ./ Y_Actual;
    
end