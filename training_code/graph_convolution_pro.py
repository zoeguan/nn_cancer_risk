from keras.layers.core import *

from keras import backend as K
from keras.engine.topology import Layer
from keras import initializers

class GraphConvPro(Layer):
    '''Subset proband's outputs from previous layer

    # Arguments
                    
    # Input shape    
        3D tensor with shape:
        `(batch_size, features, input_dim)`.
    # Output shape
        3D tensor with shape:
        `(batch_size, features, filters)`.
    '''
        
    def __init__(self,           
                 **kwargs):

        if K.backend() != 'theano':
            raise Exception("GraphConv Requires Theano Backend.")
            
        if 'input_shape' not in kwargs and 'input_dim' in kwargs:
           kwargs['input_shape'] = (kwargs.pop('input_dim'),)
           
        super(GraphConvPro, self).__init__(**kwargs)        

        self.input_spec = InputSpec(ndim=3)      

    def build(self, input_shape):
        input_dim = input_shape[2]
        self.built = True
                       
    def call(self, x):
        
        return x[:,0:1,:]

    def compute_output_shape(self, input_shape):
        return (input_shape[0], 1, input_shape[2])

