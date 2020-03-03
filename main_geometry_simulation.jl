#### Main program for sphere-tumor-geometry
include("tumor_geometry.jl")

#include("MontecarloSimulationFunctions.jl")
### Parameters that need to be defined

### Center of the tumor (x,y,z)
centerTumor = centerSphere(0,0,.15)
#centerTumor = centerSphere(0,0,.4)
### Max radius of the sphere (tumor plus spherical layers around it)
totalRad = .1
#totalRad = 343
### Number of photons
noPhotons = 100000000
### Number of layers
noLayer = 4
### Number of spherical layers
noLayerSphere = 5
### The threshold to enter the russianRoulette
thresholdWeight = .0001
### If plot 3d =1 then it will plot
p3d = 0

function main(thresholdWeight,noPhotons,noLayer,noLayerSphere,centerTumor,p3d,totalRad)
        ##Thinckness of spherical shells in microns
        microns_per_shell = 50
        ### Initialize the layers in an array
        layerArray = initializeLayersModel6()
        ### Initialize spherical layers
        layerArraySphere = initializeLayerModel6Sphere()
        ### Create arrays to plot the position of the photon when it died
        xA = zeros(noPhotons);yA = zeros(noPhotons);zA = zeros(noPhotons)
        ### cycle of each photon
        #@sync @distributed
             for i = 1:noPhotons
                ### arguments are thresholdWeight, layerArray and noLayer
                x,y,z,layerArray,layerArraySphere = simulateOnePhotonTumorGeometry(thresholdWeight,layerArray,layerArraySphere,noLayer,noLayerSphere,centerTumor,totalRad)
                ### Check if the photon is not in infinte (if it is go to zero)
                if x == -Inf || x == Inf
                    x = 0
                end
                if y == -Inf || y == Inf
                    y = 0
                end
                if z == -Inf || z == Inf
                    z = 0
                end
                xA[i] = x; yA[i] = y; zA[i] = z
            end

            if p3d == 1
                scatter(xA[:],yA[:],zA[:],title = "Position of photon when it died")
            end
            zA = sort(zA)
            return layerArray,layerArraySphere
end

layerArray,layerArraySphere = @time main(thresholdWeight,noPhotons,noLayer,noLayerSphere,centerTumor,p3d,totalRad)
energyVector = zeros(noLayer)
energyVectorSphere = zeros(noLayerSphere)

for i = 1:noLayer
    print("Layer",i," energy : ",layerArray[i].energy,"\n")
    energyVector[i] = layerArray[i].energy
end

for i = 1:noLayerSphere
    print("Layer",i," energySphere : ",layerArraySphere[i].energy,"\n")
    energyVectorSphere[i] = layerArraySphere[i].energy
end


plot()
plot!(1:noLayer,energyVector[:], title = "Energy on layers",xlabel = "Layer", ylabel = "Energy")

plot()
plot!(1:noLayerSphere,energyVectorSphere[:], title = "Energy on sphere",xlabel = "Layer", ylabel = "Energy")
energyVector = vcat(0,energyVector)
#suma = [0,0,0,0,((sum(energyVector)+sum(energyVectorSphere))*100)/noPhotons]
### add a 0 of energy in the layer of the air
data = [energyVector,energyVectorSphere]

CSV.write("DataModel_6.1.csv",DataFrame(data))
