import tensorflow as tf
import build_tvt_sets


def load_model(filepath):
    model = tf.keras.models.load_model(filepath)
    return model

#def load_model_test():
#    model1 = load_model(filepath="selu_model.h5")
#    model2 = load_model(filepath="selu_model_batch_normalized.h5")
#    val = build_tvt_sets.build_tvt()[1]
#    model1.evaluate(val[0], val[1], verbose=2)
#    model2.evaluate(val[0], val[1], verbose=2)


