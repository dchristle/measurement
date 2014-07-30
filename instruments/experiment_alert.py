# Experiment alert
# David Christle <christle@uchicago.edu> June, 2014
#
# This instrument will implement a few methods to alert via e-mail and text
# the experimenter. This could be used to indicate that a piece of equipment
# like a cryostat is outside the normal experimental parameter bounds, or just
# that a measurement is complete, for example.



from instrument import Instrument
import types
import qt
import gobject
import pprint
import time
import hdf5_data as h5
import smtplib
import numpy as np
from email.mime.text import MIMEText


class experiment_alert(Instrument):

    def __init__(self, name):
        Instrument.__init__(self, name, tags=['virtual'])
        self.name = name

        self.add_parameter('sms_email',
            flags=Instrument.FLAG_SET,
            type=types.StringType)

        self.add_parameter('email',
            flags=Instrument.FLAG_SET,
            type=types.StringType)

        self.add_parameter('gmailuser',
            flags=Instrument.FLAG_SET,
            type=types.StringType)

        self.add_parameter('gmailpw',
            flags=Instrument.FLAG_SET,
            type=types.StringType)



    def email_alert(self, message):
        server = smtplib.SMTP( "smtp.gmail.com", 587 )
        server.starttls()
        server.login( self._gmailuser, self._gmailpw )
        header = 'To:' + self._email + '\n' + 'From: ' + self._gmailuser + '\n' + 'Subject: experiment alert \n'
        msg = header + ('\n %s \n\n' % message)
        server.sendmail( self._gmailuser, self._email, msg )
        return

    def sms_alert(self, message):
        server = smtplib.SMTP( "smtp.gmail.com", 587 )
        server.starttls()
        server.login( self._gmailuser, self._gmailpw )
        header = 'To:' + self._email + '\n' + 'From: ' + self._gmailuser + '\n' + 'Subject: experiment alert \n'
        msg = header + ('\n %s \n\n' % message)
        server.sendmail( self._gmailuser, self._sms_email, msg )
        return

    def do_set_email(self, email):
        self._email = email
        return

    def do_set_sms_email(self, sms_email):
        self._sms_email = sms_email
        return
    def do_set_gmailuser(self, user):
        self._gmailuser = user
        return
    def do_set_gmailpw(self, pw):
        self._gmailpw = pw
        return
    def send_email(self,recipients,subject,message):
        '''Sends an email to 'recipients' (string or list of strings), with subject and message'''
        if isinstance(recipients,str):recipients=[recipients]
        try:
            for recipient in recipients:
                self._actually_send(recipient,subject,message)
        except smtplib.SMTPException as mailerror:
            print "Mailing error:", mailerror
            return False
        return True

    def _actually_send(self,recipient,subject,message):
        msg = MIMEText(message)
        sender = self._gmailuser
        msg['Subject'] = subject
        msg['From'] = sender
        msg['To'] = recipient

        server = smtplib.SMTP(self._host)
        server.starttls()
        server.login(self._username,self._password)
        server.sendmail(sender, recipient, msg.as_string())
        server.quit()

