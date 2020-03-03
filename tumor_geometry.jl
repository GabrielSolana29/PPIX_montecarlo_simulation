### In this file we are going to simulate the geometry proposed by Dr. Horacio
using CSV,DataFrames,StatsBase,Statistics,Distributions,Plots,Distributed
# check the radius of the photon with respect of the center (Tumor) and then decide if it is inside the
# spheres or if it somwhere else

mutable struct centerSphere
  x::Float64
  y::Float64
  z::Float64
end

### Struct is mutable so we can modify the values ###
mutable struct Photon
### Cartesian Coordinates ###
  x::Float64
  y::Float64
  z::Float64
### Directional Cosines of photon direction ###
  ux::Float64
  uy::Float64
  uz::Float64
### Weight (Fraction of the photon that goes throug layers) ###
  w::Float64
### Dead 0/1 if photon is propagating/terminated ###
  dead::Int64
### Dimensionless step size to be taken ###
### step in mm ###
  stepSize::Float64
### Number of scatterings ###
  scatters::Int64
### Index to layer wherte photon packet resides ###
  layer::Int8
end

mutable struct Layer
  ### thickness of the layer ###
  thickness::Float64
  ### refractive index ###
  refractive::Float64
  ### absorption coefficient Ma ###
  Ma::Float64
  ### scattering coefficient Ms ###
  Ms::Float64
  ### anisotropic factor g ###
  g::Float64
  ### Interaction coefficient or ###
  Mi::Float64
  ### energy stored in the layer ###
  energy::Float64
  ### x,y,zcoordinates of a layer ###
  z::Float64
  x::Float64
  y::Float64
end

mutable struct LayerSphere
  ### thickness of the layer ###
  thickness::Float64
  ### refractive index ###
  refractive::Float64
  ### absorption coefficient Ma ###
  Ma::Float64
  ### scattering coefficient Ms ###
  Ms::Float64
  ### anisotropic factor g ###
  g::Float64
  ### Interaction coefficient or ###
  Mi::Float64
  ### energy stored in the layer ###
  energy::Float64
  ### Z coordinates of a layer ###
  z::Float64
  #z1::Float32
  x::Float64
  y::Float64
end

### total radius is the radius of the total sphere
function radiusPhoton(centerSphere,photon)
  radiusPhoton = sqrt((photon.x - centerSphere.x)^2 + (photon.y - centerSphere.y)^2 + (photon.z - centerSphere.z)^2 )
  return radiusPhoton
end

function initializePhoton()
  ### Initialize x=0, y=0, z=0, ux = 0, uy = 0, uz = 1, weight = 1
  ### dead = 0 (alive), stepSize = 0, scatters = 1, layer = 0
  cartesianCoordinates = Photon(0,0,0,0,0,1,1,0,0,1,0)
  return cartesianCoordinates
end

### Decides wether a photon survives or not, given that the w <= w_th ###
### once chance in m to survive
function russianRoulette(m)
  #survive = rand(1:m)
  survive = rand(Uniform(0,1))
  if survive == m
    return true
  else
    return false
  end
end


### First step
function firstStep(layerArray)
  ### create a random number between 0-1 with a step of .001
  #x = rand(0:.0001:1)
  x = rand(Uniform(0,1))
  ### -ln(x)/ layer total interaction coefficient
  s = -log(x) /layerArray[1].Mi
  return s
end

### Generate Step size ###
function stepPhoton(currentLayer,layerArray)
  ### create a random number between 0-1 with a step of .001
  #x = rand(0:.0001:1)
  x = rand(Uniform(0,1))
  ### -ln(x)/ layer total interaction coefficient
  s = -log(x) / (layerArray[currentLayer].Mi)
  return s
end

### Move photon
function move(photon)
  ### Update positions
  photon.x = photon.x + (photon.ux * photon.stepSize)
  photon.y = photon.y + (photon.uy * photon.stepSize)
  photon.z = photon.z + (photon.uz * photon.stepSize)
  return photon
end


function spin(currentLayer,layerArray,photon)
  ### Sample for costheta
  #rnd = rand(0:.0001:1)
  rnd = rand(Uniform(0,1))
  g = layerArray[currentLayer].g
  costheta = 0
  if g == 0
    costheta = (2*rnd) - 1
  else
    temp = (1-(g*g)) / (1- g + (2*g*rnd))
    costheta = (1 + (g*g) - (temp*temp)) / (2*g)
  end
################################## Problema ##############################
#####   En las simulaciones me esta dando un numero negativo por lo que la raiz se convierte en compleja    #####
  #print(" " ,costheta)
  ayuda = 1 - (costheta * costheta)
  if ayuda < .000000000001
    ayuda = 0
  end
  sintheta = sqrt(ayuda) ## sqrt is faster than sin()
  ### Sample psi
  rnd = rand(0:.0001:1)
  psi = 2*pi*rnd
  cospsi = cos(psi) ### is radians

  if psi < pi
    sinpsi = sqrt(1-(cospsi*cospsi))
  else
    sinpsi = -(sqrt(1-(cospsi*cospsi)))
  end
  #### New trajectory
  if 1 - abs(photon.uz) < .000000000001 ###
    uxx = sintheta * cospsi
    uyy = sintheta * sinpsi
    uzz = costheta * sign(photon.uz) ### sign returns the sign of the number (+1 or -1)
  else
    temp = sqrt(1 - (photon.uz * photon.uz))
    uxx = (sintheta * ((photon.ux * photon.uz * cospsi) - (photon.uy * sinpsi)) / temp) + (photon.ux * costheta)
    uyy = (sintheta * ((photon.uy * photon.uz * cospsi) + (photon.ux * sinpsi)) / temp) + (photon.uy * costheta)
    uzz = (-sintheta * cospsi * temp) + (photon.uz * costheta)
  end
  ### update trajectory
  photon.ux = uxx
  photon.uy = uyy
  photon.uz = uzz
  return photon
end

function absorb(photon,currentLayer,layerArray)
    ### This parameter depends on the layer
    ### Albedo equals the fractional probability of being scattered
    albedo = layerArray[currentLayer].Ms / layerArray[currentLayer].Mi
    ### Photon weight absorbed at this step
    absorbt = photon.w * (1-albedo)
    #absorbt = .9817 * (1-albedo)
    ### Decrement weight  by amount absorbed
    photon.w = photon.w - absorbt
    ### Save the energy that is placed in the layer
    layerArray[currentLayer].energy = layerArray[currentLayer].energy + absorbt
    ### update bins
    ####Falta implementar, platicar con Dr. horacio
    ### update the different bins (spherical, cylindrical, and planar arryas)
    return photon
end

### Return the current layer where the particle is dropping energy
function currentLayerUpdate(layerArray,noLayer,photon)
    distZ = 0
    ### Through the cycle you check in which layer the photon.z is
    for i = 1:noLayer
      distZ = distZ + layerArray[i].thickness
      if photon.z <= distZ && photon.z > 0
        return i
      elseif photon.z <= 0
        return 0
      ##elseif photon.z > layerArray[i].thickness && i == noLayer
    elseif photon.z > distZ && i == noLayer
        return 0
      end
    end
    return 0
end

function currentLayerUpdateSpherical(layerArraySphere,noLayerSphere,radPhoton,totalRad)
  radLayer = totalRad
  radLayer_2 = 0
  for i = 1:noLayerSphere
    if i == noLayerSphere
      return noLayerSphere
    end
    radLayer = radLayer - (radLayer - (layerArraySphere[i].thickness/2))
    radLayer_2 = radLayer - (radLayer - (layerArraySphere[i+1].thickness/2))
    if radPhoton <= radLayer && radPhoton > radLayer_2
      return i
    end
  end
end

### Initilize simulation
function simulateOnePhotonTumorGeometry(thersholdWeight,layerArray,layerArraySphere,noLayer,noLayerSphere,centerTumor,totalRad)
  currentLayer = 1
  currentLayerSphere = 1
  ### Create a new photon
  photon = initializePhoton() #se queda igual
  sphereLinear = 0
  ### While photon is not dead
  while photon.dead == 0
    ### if there is no stepSize ###
    if photon.stepSize == 0
      #costheta = 2*rand(0:.0001:1) - 1
      costheta = (2*rand(Uniform(0,1))) - 1
      sintheta = sqrt(1-(costheta*costheta))
      #psi = 2*pi*rand(0:.0001:1)
      psi = 2*pi*rand(Uniform(0,1))
      photon.ux = sintheta * cos(psi)
      photon.uy = sintheta * sin(psi)
      photon.uz = costheta
      photon.stepSize = firstStep(layerArray) #se queda igual
    elseif sphereLinear == 0
      photon.stepSize = stepPhoton(currentLayer,layerArray) #se queda igual
    elseif sphereLinear == 1
      photon.stepSize = stepPhoton(currentLayerSphere,layerArraySphere) #se queda igual
    end
    sphereLinear = 0
    ### update the photon position
    photon = move(photon) #se queda igual
    ### Check if the photon is now in a new layer, inside the sphere or not)
    radPhoton = radiusPhoton(centerTumor,photon)
    ### If radPhoton is bigger than total Rad it means is not in the sphere
    if radPhoton <= totalRad
      sphereLinear = 1
      currentLayerSphere = currentLayerUpdateSpherical(layerArraySphere,noLayerSphere,radPhoton,totalRad)
      ### Update the weight of the photon after absorption of the layer
      photon = absorb(photon,currentLayerSphere,layerArraySphere)
      ### update the current trajectory ux,uy,uz directional cosines
      photon = spin(currentLayerSphere,layerArraySphere,photon) # completo
    elseif radPhoton > totalRad && abs(photon.x) < layerArray[currentLayer].x && abs(photon.y) < layerArray[currentLayer].y
      currentLayer = currentLayerUpdate(layerArray,noLayer,photon)
      ### If currentLayer is == 0 means that the photon was reflected and went out of the layers
      ### So we kill it
      if currentLayer == 0
        return photon.x,photon.y,photon.z,layerArray,layerArraySphere
      end
      ### Update the weight of the photon after absorption of the layer
      photon = absorb(photon,currentLayer,layerArray) # falta bins
      ### update the current trajectory ux,uy,uz directional cosines
      photon = spin(currentLayer,layerArray,photon) # completo
    else
      return photon.x,photon.y,photon.z,layerArray,layerArraySphere
    end

    #scatter()
    ### russian roulette gives 1/10 chances to survive after the photon has
    ### .00001 weight left
    if photon.w < thresholdWeight
      if russianRoulette(10) == true
        photon.w == photon.w /.1
      else
        photon.w = 0
        photon.dead = 1
        return photon.x,photon.y,photon.z,layerArray,layerArraySphere
      end
    end
  end
end


#### Modelo 1
function initializeLayersModel1()
  ### Initialize layer struct for the first layer
  ### parameters matching the tiny.c code to compare
  ### thickness=.1,refractive=.1,Ma=.1,Ms=.1,g=.1,Mi=0,energy=0,z0=0,x = ,y =###
  ### air
  #layer1 = Layer(.8,1,0,1,1,0,0,0,20,20)
  #layer1.Mi = layer1.Ma + layer1.Ms

  ### layer skin 1
  layer2 = Layer(.3,1.38,.7,36.7,.852,0,0,0,1,1)
  layer2.Mi = layer2.Ma + layer2.Ms

  ### layer lineal skin
  layer3 = Layer(.3,1.38,.7,36.7,.852,0,0,0,1,1)
  layer3.Mi = layer3.Ma + layer3.Ms

  ### layer lineal skin
  layer4 = Layer(.3,1.38,.7,36.7,.852,0,0,0,1,1)
  layer4.Mi = layer4.Ma + layer4.Ms

  ### layer lineal skin
  layer5 = Layer(.3,1.38,.7,36.7,.852,0,0,0,1,1)
  layer5.Mi = layer5.Ma + layer5.Ms
  ### Create an array of struct Layer of 1 dimension
  layerArray = Array{Layer,1}(undef,0)

  ### Push layers in to array
  #push!(layerArray,layer1)
  push!(layerArray,layer2)
  push!(layerArray,layer3)
  push!(layerArray,layer4)
  push!(layerArray,layer5)
  return layerArray
end


function initializeLayerModel1Sphere()
  ### the spheres are one inside the other
   ### layer skin spherical layer 4
   layer1 = LayerSphere(.2,1.38,.7,36.7,.852,0,0,0,0,0)
   layer1.Mi = layer1.Ma + layer1.Ms

   ### layer  skin spherical layer 3
   layer2 = LayerSphere(.175,1.38,.7,36.7,.852,0,0,0,0,0)
   layer2.Mi = layer2.Ma + layer2.Ms

   ### layer skin spherical layer 2
   layer3 = LayerSphere(.15,1.38,.7,36.7,.852,0,0,0,0,0)
   layer3.Mi = layer3.Ma + layer3.Ms

   ### layer skin spherical layer1
   layer4 = LayerSphere(.125,1.38,.7,36.7,.852,0,0,0,0,0)
   layer4.Mi = layer4.Ma + layer4.Ms

   ### layer spherical tumor
   layer5 = LayerSphere(.1,1.38,2.3,21.2,.852,0,0,0,0,0)
   layer5.Mi = layer5.Ma + layer5.Ms

   ### Create an array of struct Layer of 1 dimension
   layerArraySphere = Array{LayerSphere,1}(undef,0)
   ### Push layers in to array
   push!(layerArraySphere,layer1)
   push!(layerArraySphere,layer2)
   push!(layerArraySphere,layer3)
   push!(layerArraySphere,layer4)
   push!(layerArraySphere,layer5)
   return layerArraySphere
end

#### Modelo 2
function initializeLayersModel2()
  ### Initialize layer struct for the first layer
  ### parameters matching the tiny.c code to compare
  ### thickness=.1,refractive=.1,Ma=.1,Ms=.1,g=.1,Mi=0,energy=0,z0=0,x = ,y =###
  ### air
  #layer1 = Layer(.8,1,0,1,1,0,0,0,20,20)
  #layer1.Mi = layer1.Ma + layer1.Ms

  ### layer skin 1
  layer2 = Layer(.3,1.38,.7,36.7,.852,0,0,0,1,1)
  layer2.Mi = layer2.Ma + layer2.Ms

  ### layer lineal skin
  layer3 = Layer(.3,1.38,.7,36.7,.852,0,0,0,1,1)
  layer3.Mi = layer3.Ma + layer3.Ms

  ### layer lineal skin
  layer4 = Layer(.3,1.38,.7,36.7,.852,0,0,0,1,1)
  layer4.Mi = layer4.Ma + layer4.Ms

  ### layer lineal skin
  layer5 = Layer(.3,1.38,.7,36.7,.852,0,0,0,1,1)
  layer5.Mi = layer5.Ma + layer5.Ms
  ### Create an array of struct Layer of 1 dimension
  layerArray = Array{Layer,1}(undef,0)

  ### Push layers in to array
  #push!(layerArray,layer1)
  push!(layerArray,layer2)
  push!(layerArray,layer3)
  push!(layerArray,layer4)
  push!(layerArray,layer5)
  return layerArray
end


function initializeLayerModel2Sphere()
  ### the spheres are one inside the other
   ### layer skin spherical layer 4
   layer1 = LayerSphere(.2,1.38,.7,36.7,.852,0,0,0,0,0)
   layer1.Mi = layer1.Ma + layer1.Ms

   ### layer  skin spherical layer 3
   layer2 = LayerSphere(.175,1.38,.7,36.7,.852,0,0,0,0,0)
   layer2.Mi = layer2.Ma + layer2.Ms

   ### layer skin spherical layer 2
   layer3 = LayerSphere(.15,1.38,.7,36.7,.852,0,0,0,0,0)
   layer3.Mi = layer3.Ma + layer3.Ms

   ### layer skin spherical layer1
   layer4 = LayerSphere(.125,1.38,.7,36.7,.852,0,0,0,0,0)
   layer4.Mi = layer4.Ma + layer4.Ms

   ### layer spherical tumor
   layer5 = LayerSphere(.1,1.38,14.3,22.4,.852,0,0,0,0,0)
   layer5.Mi = layer5.Ma + layer5.Ms

   ### Create an array of struct Layer of 1 dimension
   layerArraySphere = Array{LayerSphere,1}(undef,0)
   ### Push layers in to array
   push!(layerArraySphere,layer1)
   push!(layerArraySphere,layer2)
   push!(layerArraySphere,layer3)
   push!(layerArraySphere,layer4)
   push!(layerArraySphere,layer5)
   return layerArraySphere
end

#### Modelo 3
function initializeLayersModel3()
  ### Initialize layer struct for the first layer
  ### parameters matching the tiny.c code to compare
  ### thickness=.1,refractive=.1,Ma=.1,Ms=.1,g=.1,Mi=0,energy=0,z0=0,x = ,y =###
  ### air
  #layer1 = Layer(.8,1,0,1,1,0,0,0,20,20)
  #layer1.Mi = layer1.Ma + layer1.Ms

  ### layer skin 1
  layer2 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer2.Mi = layer2.Ma + layer2.Ms

  ### layer lineal skin
  layer3 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer3.Mi = layer3.Ma + layer3.Ms

  ### layer lineal skin
  layer4 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer4.Mi = layer4.Ma + layer4.Ms

  ### layer lineal skin
  layer5 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer5.Mi = layer5.Ma + layer5.Ms
  ### Create an array of struct Layer of 1 dimension
  layerArray = Array{Layer,1}(undef,0)

  ### Push layers in to array
  #push!(layerArray,layer1)
  push!(layerArray,layer2)
  push!(layerArray,layer3)
  push!(layerArray,layer4)
  push!(layerArray,layer5)
  return layerArray
end


function initializeLayerModel3Sphere()
  ### the spheres are one inside the other
   ### layer skin spherical layer 4
   layer1 = LayerSphere(.2,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer1.Mi = layer1.Ma + layer1.Ms

   ### layer  skin spherical layer 3
   layer2 = LayerSphere(.175,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer2.Mi = layer2.Ma + layer2.Ms

   ### layer skin spherical layer 2
   layer3 = LayerSphere(.15,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer3.Mi = layer3.Ma + layer3.Ms

   ### layer skin spherical layer1
   layer4 = LayerSphere(.125,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer4.Mi = layer4.Ma + layer4.Ms

   ### layer spherical tumor
   layer5 = LayerSphere(.1,1.38,2.3,21.2,.8027,0,0,0,0,0)
   layer5.Mi = layer5.Ma + layer5.Ms

   ### Create an array of struct Layer of 1 dimension
   layerArraySphere = Array{LayerSphere,1}(undef,0)
   ### Push layers in to array
   push!(layerArraySphere,layer1)
   push!(layerArraySphere,layer2)
   push!(layerArraySphere,layer3)
   push!(layerArraySphere,layer4)
   push!(layerArraySphere,layer5)
   return layerArraySphere
end

#### Modelo 2
function initializeLayersModel4()
  ### Initialize layer struct for the first layer
  ### parameters matching the tiny.c code to compare
  ### thickness=.1,refractive=.1,Ma=.1,Ms=.1,g=.1,Mi=0,energy=0,z0=0,x = ,y =###
  ### air
  #layer1 = Layer(.8,1,0,1,1,0,0,0,20,20)
  #layer1.Mi = layer1.Ma + layer1.Ms

  ### layer skin 1
  layer2 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer2.Mi = layer2.Ma + layer2.Ms

  ### layer lineal skin
  layer3 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer3.Mi = layer3.Ma + layer3.Ms

  ### layer lineal skin
  layer4 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer4.Mi = layer4.Ma + layer4.Ms

  ### layer lineal skin
  layer5 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer5.Mi = layer5.Ma + layer5.Ms
  ### Create an array of struct Layer of 1 dimension
  layerArray = Array{Layer,1}(undef,0)

  ### Push layers in to array
  #push!(layerArray,layer1)
  push!(layerArray,layer2)
  push!(layerArray,layer3)
  push!(layerArray,layer4)
  push!(layerArray,layer5)
  return layerArray
end


function initializeLayerModel4Sphere()
  ### the spheres are one inside the other
   ### layer skin spherical layer 4
   layer1 = LayerSphere(.2,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer1.Mi = layer1.Ma + layer1.Ms

   ### layer  skin spherical layer 3
   layer2 = LayerSphere(.175,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer2.Mi = layer2.Ma + layer2.Ms

   ### layer skin spherical layer 2
   layer3 = LayerSphere(.15,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer3.Mi = layer3.Ma + layer3.Ms

   ### layer skin spherical layer1
   layer4 = LayerSphere(.125,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer4.Mi = layer4.Ma + layer4.Ms

   ### layer spherical tumor
   layer5 = LayerSphere(.1,1.5,79.28,28.8981,.8027,0,0,0,0,0)
   layer5.Mi = layer5.Ma + layer5.Ms

   ### Create an array of struct Layer of 1 dimension
   layerArraySphere = Array{LayerSphere,1}(undef,0)
   ### Push layers in to array
   push!(layerArraySphere,layer1)
   push!(layerArraySphere,layer2)
   push!(layerArraySphere,layer3)
   push!(layerArraySphere,layer4)
   push!(layerArraySphere,layer5)
   return layerArraySphere
end


#### Modelo 4.2 with small concentration of PPIX 69microMolar
function initializeLayersModel5()
  ### Initialize layer struct for the first layer
  ### parameters matching the tiny.c code to compare
  ### thickness=.1,refractive=.1,Ma=.1,Ms=.1,g=.1,Mi=0,energy=0,z0=0,x = ,y =###
  ### air
  #layer1 = Layer(.8,1,0,1,1,0,0,0,20,20)
  #layer1.Mi = layer1.Ma + layer1.Ms

  ### layer skin 1
  layer2 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer2.Mi = layer2.Ma + layer2.Ms

  ### layer lineal skin
  layer3 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer3.Mi = layer3.Ma + layer3.Ms

  ### layer lineal skin
  layer4 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer4.Mi = layer4.Ma + layer4.Ms

  ### layer lineal skin
  layer5 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer5.Mi = layer5.Ma + layer5.Ms
  ### Create an array of struct Layer of 1 dimension
  layerArray = Array{Layer,1}(undef,0)

  ### Push layers in to array
  #push!(layerArray,layer1)
  push!(layerArray,layer2)
  push!(layerArray,layer3)
  push!(layerArray,layer4)
  push!(layerArray,layer5)
  return layerArray
end


function initializeLayerModel5Sphere()
  ### the spheres are one inside the other
   ### layer skin spherical layer 4
   layer1 = LayerSphere(.2,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer1.Mi = layer1.Ma + layer1.Ms

   ### layer  skin spherical layer 3
   layer2 = LayerSphere(.175,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer2.Mi = layer2.Ma + layer2.Ms

   ### layer skin spherical layer 2
   layer3 = LayerSphere(.15,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer3.Mi = layer3.Ma + layer3.Ms

   ### layer skin spherical layer1
   layer4 = LayerSphere(.125,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer4.Mi = layer4.Ma + layer4.Ms

   ### layer spherical tumor
   layer5 = LayerSphere(.1,1.5,3.1136,21.2813,.8027,0,0,0,0,0)
   layer5.Mi = layer5.Ma + layer5.Ms

   ### Create an array of struct Layer of 1 dimension
   layerArraySphere = Array{LayerSphere,1}(undef,0)
   ### Push layers in to array
   push!(layerArraySphere,layer1)
   push!(layerArraySphere,layer2)
   push!(layerArraySphere,layer3)
   push!(layerArraySphere,layer4)
   push!(layerArraySphere,layer5)
   return layerArraySphere
end


#### Modelo 4.2 with small concentration of PPIX 69microMolar
function initializeLayersModel6()
  ### Initialize layer struct for the first layer
  ### parameters matching the tiny.c code to compare
  ### thickness=.1,refractive=.1,Ma=.1,Ms=.1,g=.1,Mi=0,energy=0,z0=0,x = ,y =###
  ### air
  #layer1 = Layer(.8,1,0,1,1,0,0,0,20,20)
  #layer1.Mi = layer1.Ma + layer1.Ms

  ### layer skin 1
  layer2 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer2.Mi = layer2.Ma + layer2.Ms

  ### layer lineal skin
  layer3 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer3.Mi = layer3.Ma + layer3.Ms

  ### layer lineal skin
  layer4 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer4.Mi = layer4.Ma + layer4.Ms

  ### layer lineal skin
  layer5 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer5.Mi = layer5.Ma + layer5.Ms
  ### Create an array of struct Layer of 1 dimension
  layerArray = Array{Layer,1}(undef,0)

  ### Push layers in to array
  #push!(layerArray,layer1)
  push!(layerArray,layer2)
  push!(layerArray,layer3)
  push!(layerArray,layer4)
  push!(layerArray,layer5)
  return layerArray
end


function initializeLayerModel6Sphere()
  ### the spheres are one inside the other
   ### layer skin spherical layer 4
   layer1 = LayerSphere(.2,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer1.Mi = layer1.Ma + layer1.Ms

   ### layer  skin spherical layer 3
   layer2 = LayerSphere(.175,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer2.Mi = layer2.Ma + layer2.Ms

   ### layer skin spherical layer 2
   layer3 = LayerSphere(.15,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer3.Mi = layer3.Ma + layer3.Ms

   ### layer skin spherical layer1
   layer4 = LayerSphere(.125,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer4.Mi = layer4.Ma + layer4.Ms

   ### layer spherical tumor
   layer5 = LayerSphere(.1,1.5,3.4791,21.32,.8027,0,0,0,0,0)
   layer5.Mi = layer5.Ma + layer5.Ms

   ### Create an array of struct Layer of 1 dimension
   layerArraySphere = Array{LayerSphere,1}(undef,0)
   ### Push layers in to array
   push!(layerArraySphere,layer1)
   push!(layerArraySphere,layer2)
   push!(layerArraySphere,layer3)
   push!(layerArraySphere,layer4)
   push!(layerArraySphere,layer5)
   return layerArraySphere
end


#### Modelo 4.2 with small concentration of PPIX 69microMolar
function initializeLayersModel7()
  ### Initialize layer struct for the first layer
  ### parameters matching the tiny.c code to compare
  ### thickness=.1,refractive=.1,Ma=.1,Ms=.1,g=.1,Mi=0,energy=0,z0=0,x = ,y =###
  ### air
  #layer1 = Layer(.8,1,0,1,1,0,0,0,20,20)
  #layer1.Mi = layer1.Ma + layer1.Ms

  ### layer skin 1
  layer2 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer2.Mi = layer2.Ma + layer2.Ms

  ### layer lineal skin
  layer3 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer3.Mi = layer3.Ma + layer3.Ms

  ### layer lineal skin
  layer4 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer4.Mi = layer4.Ma + layer4.Ms

  ### layer lineal skin
  layer5 = Layer(.3,1.38,.7,36.7,.8027,0,0,0,1,1)
  layer5.Mi = layer5.Ma + layer5.Ms
  ### Create an array of struct Layer of 1 dimension
  layerArray = Array{Layer,1}(undef,0)

  ### Push layers in to array
  #push!(layerArray,layer1)
  push!(layerArray,layer2)
  push!(layerArray,layer3)
  push!(layerArray,layer4)
  push!(layerArray,layer5)
  return layerArray
end


function initializeLayerModel7Sphere()
  ### the spheres are one inside the other
   ### layer skin spherical layer 4
   layer1 = LayerSphere(.2,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer1.Mi = layer1.Ma + layer1.Ms

   ### layer  skin spherical layer 3
   layer2 = LayerSphere(.175,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer2.Mi = layer2.Ma + layer2.Ms

   ### layer skin spherical layer 2
   layer3 = LayerSphere(.15,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer3.Mi = layer3.Ma + layer3.Ms

   ### layer skin spherical layer1
   layer4 = LayerSphere(.125,1.38,.7,36.7,.8027,0,0,0,0,0)
   layer4.Mi = layer4.Ma + layer4.Ms

   ### layer spherical tumor
   layer5 = LayerSphere(.1,1.5,14.092,22.38,.8027,0,0,0,0,0)
   layer5.Mi = layer5.Ma + layer5.Ms

   ### Create an array of struct Layer of 1 dimension
   layerArraySphere = Array{LayerSphere,1}(undef,0)
   ### Push layers in to array
   push!(layerArraySphere,layer1)
   push!(layerArraySphere,layer2)
   push!(layerArraySphere,layer3)
   push!(layerArraySphere,layer4)
   push!(layerArraySphere,layer5)
   return layerArraySphere
end
