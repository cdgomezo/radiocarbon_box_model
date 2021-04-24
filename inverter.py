#Author: Carlos Gómez-Ortiz
#Affiliation: Department of Physical Geography and Ecosystem Science, Lund University
#2021

from numpy import *
from numpy.linalg import *
from copy import *

class Inverter:
    def __init__(self,
                 trueFluxes,
                 trueInit,
                 model,
                 presFluxKey = [],
                 realObs = False,
                 fromTrue = {},
                 trueObs = None,
                 realError = False,
                 errorObs = {},
                 errorInit = {},
                 errorFluxes = {}):
        """
        Class that contains the functions necessary to convert from flux to vector and viceversa,
        observations to vector and viceversa, \Delta^14CO_2 to F^14CO2 and viceversa, and the inversion:
        Kernel matrix, Gain matrix, Cost function, and the posterior computation.
        """
        
        self.trueFluxes = trueFluxes
        self.trueInit = trueInit
        self.model = model
        self.presFluxKey = presFluxKey
        self.realObs = realObs
        self.fromTrue = fromTrue
        self.trueObs = trueObs
        self.realError = realError
        self.errorObs = errorObs
        self.errorInit = errorInit
        self.errorFluxes = errorFluxes
       
        
    def genPrior(self): 
        """
        Generates prior fluxes from true fluxes and returns the control vector.
        """
        self.priorFluxes = deepcopy(self.trueFluxes)
        self.trueVector = self.flux2Vector(self.trueInit, self.trueFluxes)
        
        if self.realObs:
            self.controlVector = self.flux2Vector(self.trueInit, self.priorFluxes)
        else:
            for key_i, value_i in self.priorFluxes.items():
                if key_i == 'time' or key_i == 'units':
                    pass
                else:
                    for key_j, value_j in self.priorFluxes[key_i].items():
                        for key_k, value_k in self.priorFluxes[key_i][key_j].items():
                            for key_l, value_l in self.presFluxKey.items():
                                if key_i == key_l:
                                    if key_k in self.presFluxKey[key_l]:
                                        pass
                                    else:
                                        self.priorFluxes[key_i][key_j][key_k] = value_k * self.fromTrue[key_i][key_j][key_k]

            self.controlVector = self.flux2Vector(self.trueInit, self.priorFluxes) 

        return self.priorFluxes, self.controlVector, self.trueVector
    
    
    def genPrescribed(self):
        """
        Generates prescribed fluxes. The initial condition is prescribed by default.
        """
        prescribedInit = deepcopy(self.trueInit)

        for key_i, value_i in prescribedInit.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in prescribedInit[key_i].items():
                    prescribedInit[key_i][key_j] = value_j * 0

        self.prescribedFlux = deepcopy(self.trueFluxes)

        for key_i, value_i in self.prescribedFlux.items():
            if key_i == 'time' or key_i == 'units':
                pass
            else:
                for key_j, value_j in self.prescribedFlux[key_i].items():
                    for key_k, value_k in self.prescribedFlux[key_i][key_j].items():
                        for key_l, value_l in self.presFluxKey.items():
                            if key_i == key_l:
                                if key_k in self.presFluxKey[key_l]:
                                    pass
                                else:
                                    self.prescribedFlux[key_i][key_j][key_k] = value_k * 0

        if self.realObs:
            self.prescribedConc = self.model(prescribedInit, self.prescribedFlux)
            self.prescribedConc = self.extractObs(self.prescribedConc)
        else:
            self.prescribedConc = self.model(prescribedInit, self.prescribedFlux)
            self.prescribedConc = self.extractObs(self.prescribedConc)

        self.prescribedConcVector = self.obs2Vector(self.prescribedConc)

        return self.prescribedFlux, self.prescribedConcVector

    
    def flux2Vector(self, newInit, fluxes): 
        """
        Returns a vector with the structure of the control vector.
        """
        vector = array(())

        for key_i, value_i in newInit.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in newInit[key_i].items():
                    vector = append(vector, value_j)

        for key_i, value_i in fluxes.items():
            if key_i == 'time' or key_i == 'units':
                pass
            else:
                for key_j, value_j in fluxes[key_i].items():
                    for key_k, value_k in fluxes[key_i][key_j].items():
                        for key_l, value_l in self.presFluxKey.items():
                            if key_i == key_l:
                                if key_k in self.presFluxKey[key_l]:
                                    pass
                                else:
                                    vector = append(vector, value_k)
        
        return vector


    def vector2Flux(self, vector, prescribed = True): 
        """
        Returns flux and initial condition from a vector.
        """        
        newInit = deepcopy(self.trueInit)

        fluxes = deepcopy(self.trueFluxes)

        i = 0
        for key_i, value_i in newInit.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in newInit[key_i].items():
                    newInit[key_i][key_j] = vector[i]
                    i += 1

        j = 0
        for key_i, value_i in fluxes.items():
            if key_i == 'time' or key_i == 'units':
                pass
            else:
                for key_j, value_j in fluxes[key_i].items():
                    for key_k, value_k in fluxes[key_i][key_j].items():
                        for key_l, value_l in self.presFluxKey.items():
                            if key_i == key_l:
                                if key_k in self.presFluxKey[key_l]:
                                    if prescribed:
                                        fluxes[key_i][key_j][key_k] = value_k * 0.
                                    else:
                                        pass
                                else:
                                    fluxes[key_i][key_j][key_k] = vector[i+j*(len(fluxes['time'])):i+(j+1)*(len(fluxes['time']))].reshape(-1,1)
                                    j += 1      
        
        return newInit, fluxes


    def genObsVector(self): 
        """
        Returns the true observation vector.
        """
        if self.realObs:
            try:
                self.obsVector = self.obs2Vector(self.trueObs)
            except ValueError:
                if self.trueObs == None:
                    print("An observation set must be enter when 'realObs' is set to True.")
                else:
                    print("Observation set with invalid format.")
        else:
            self.trueObs = self.model(self.trueInit, self.trueFluxes)
            self.obsVector = self.obs2Vector(self.trueObs)

        return self.trueObs, self.obsVector


    def obs2Vector(self, obs):
        """
        Returns a vector with the structure of the observation vector.
        """
        obsVector = array(())

        for key_i, value_i in obs.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in obs[key_i].items():
                    for key_k in obs[key_i][key_j].keys():
                        if key_k == 'time':
                            pass
                        else:
                            obsVector = append(obsVector, obs[key_i][key_j][key_k])

        return obsVector


    def vector2Obs(self, vector):
        """
        Returns observation vector as observations. It is never used in the code, just for comparison with obs2Vector.
        """
        obs = deepcopy(self.trueObs)

        i = 0
        for key_i, value_i in obs.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in obs[key_i].items():
                    obs[key_i][key_j]['value'] = vector[i * len(obs[key_i][key_j]['time']):(i+1) * len(obs[key_i][key_j]['time'])]
                    i += 1

        return obs
    
    
    def genErrorVector(self):
        """
        Generates control and observation error vectors.
        """
        if self.realError:
            self.errorControlVector = self.flux2Vector(self.errorInit, self.errorFluxes)
            self.errorObsVector = self.obs2Vector(self.errorObs)

        else:
            errInit = deepcopy(self.errorInit)
            errFlux = deepcopy(self.errorFluxes)
            errObs  = deepcopy(self.errorObs)
            
            for key_i, value_i in self.trueInit.items():
                if key_i == 'units':
                    pass
                else:
                    for key_j, value_j in self.trueInit[key_i].items():
                        errInit[key_i][key_j] = value_j * self.errorInit[key_i][key_j]

            for key_i, value_i in self.priorFluxes.items():
                if key_i == 'time' or key_i == 'units':
                    pass
                else:
                    for key_j, value_j in self.priorFluxes[key_i].items():
                        for key_k, value_k in self.priorFluxes[key_i][key_j].items():
                                errFlux[key_i][key_j][key_k] = value_k * self.errorFluxes[key_i][key_j][key_k]

            self.errorControlVector = self.flux2Vector(errInit, errFlux)

            for key_i, value_i in self.trueObs.items():
                if key_i == 'units':
                    pass
                else:
                    for key_j, value_j in self.trueObs[key_i].items():
                        for key_k, value_k in self.trueObs[key_i][key_j].items():
                            if key_k == 'time':
                                pass
                            else:
                                errObs[key_i][key_j][key_k] = value_k * self.errorObs[key_i][key_j][key_k]

            self.errorObsVector = self.obs2Vector(errObs)

        return self.errorControlVector, self.errorObsVector


    def extractObs(self, obs):
        """
        Extracts the observations generated by the transport model when the number of observations is inferior than the number of fluxes.
        """

        for key_i, value_i in obs.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in obs[key_i].items():
                    o = array(())
                    obs[key_i][key_j]['time'] = self.trueObs[key_i][key_j]['time']
                    for i in range(len(self.trueFluxes['time'])):
                        if self.trueFluxes['time'][i] in self.trueObs[key_i][key_j]['time']:
                            o = append(o, obs[key_i][key_j]['value'][i])
                            if i == len(self.trueFluxes['time'])-1:
                                obs[key_i][key_j]['value'] = o.reshape(-1,1)
                        else:
                            pass

        return obs
        
        
    def cost(self, posteriorVector):
        """
        Returns the cost function.
        """
        dx = posteriorVector - self.controlVector
        dy = matmul(self.H, posteriorVector) - (self.obsVector - self.prescribedConcVector)
        
        return matmul(matmul(dx.T, inv(self.B)), dx) * 0.5 + matmul(matmul(dy.T, inv(self.R)), dy) * 0.5

    
    def inversion(self):
        """
        Contains the construction of H and G matrices, error covariance matrices B and R,
        and returns the values of the posterior, cost function and matrices.
        """
        self.B = eye(len(self.controlVector))*(self.errorControlVector**2) # Error covariance matrix for control vector
        
        self.R = eye(len(self.obsVector))*(self.errorObsVector**2) # Error covariance matrix for obs
        
        self.H = zeros((len(self.obsVector), len(self.controlVector))) # Observation operator
        for i in range(len(self.controlVector)):
            vec = zeros((len(self.controlVector)))
            vec[i] = 1
            c01, flux1 = self.vector2Flux(vec)
            conc1 = self.model(c01, flux1)
            if self.realObs:
                conc1 = self.extractObs(conc1)
            cvec = self.obs2Vector(conc1)
            self.H[:, i] = cvec

        G = matmul(matmul(self.B, self.H.T), inv(matmul(matmul(self.H, self.B), self.H.T) + self.R))
        
        posteriorVector = self.controlVector + matmul(G, ((self.obsVector - self.prescribedConcVector) - matmul(self.H, self.controlVector)))
        posteriorInit, self.posteriorFluxes = self.vector2Flux(posteriorVector, False)
        concPosterior = self.model(posteriorInit, self.posteriorFluxes)
        if self.realObs:
            concPosterior = self.extractObs(concPosterior)
        
        # Cost function values
        costTrue = self.cost(self.trueVector)
        costPrior = self.cost(self.controlVector)
        costPosterior = self.cost(posteriorVector)
        
        if self.realObs:
            print('Cost function prior:', costPrior, 'Cost function posterior:', costPosterior)
        else:
            print('Cost function true:', costTrue,'Cost function prior:', costPrior, 'Cost function posterior:', costPosterior)

        concPrior = self.model(self.trueInit, self.priorFluxes)
        if self.realObs:
            concPrior = self.extractObs(concPrior)
        
        return posteriorVector, self.posteriorFluxes, concPrior, concPosterior, self.B, self.R, self.H, G
