from lib.network.objectsharernew import objectsharer as objshn


objshn.helper.backend.start_server('127.0.0.1')
objshn.helper.backend.connect_to('tcp://205.208.82.175:12050')

st0_h = objshn.find_object('standa0')