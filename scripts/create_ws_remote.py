from lib.network import object_sharer as objsh

if objsh.start_glibtcp_client('205.208.82.175',port=12002, nretry=3):
    remote_ins_server=objsh.helper.find_object('instrument_server')
    ws = qt.instruments.create('ws_here', 'Remote_Instrument',
                 remote_name='ws', inssrv=remote_ins_server)
