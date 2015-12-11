""" Notification helpers """

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.utils import COMMASPACE, formatdate
import os
import smtplib

SMS_ADDRS = {
        "TELUS": "{0}@msg.telus.com",
        "BELL": "{0}@txt.bellmobility.ca",
        "ROGERS": "{0}@pcs.rogers.com",
        "VERIZON": "{0}@vtext.com",
        "ATT": "{0}@txt.att.net",
        "T-MOBILE": "{0}@tmomail.com",
        "SPRINT": "{0}@messaging.sprintpcs.com",
        }

def send_sms(mailconf, from_addr, number, provider, msg, subject=""):
    if provider not in SMS_ADDRS:
        raise KeyError("provider must be one of:", list(SMS_ADDRS.keys()))
    raise NotImplementedError


def send_mail(mailconf, from_addr, to_addr, msg, subject="", attachment=None):
    """ reads login information from a text file containing:

        smtp.server
        port
        username
        password

    and sends an email containing *msg*
    """

    if not os.path.isfile(mailconf):
        raise IOError("no file %s containing mail login info" % mailconf)

    with open(mailconf, "r") as f:
        host = f.readline().strip()
        port = f.readline().strip()
        user = f.readline().strip()
        passwd = f.readline().strip()

    print("%s:%s" % (host, port))
    print(user)
    print("--------------------")

    mime_msg = MIMEMultipart()
    mime_msg["From"] = from_addr
    mime_msg["To"] = to_addr
    mime_msg["Subject"] = subject
    mime_msg["Date"] = formatdate(localtime=True)
    mime_msg.attach(MIMEText(msg))

    if attachment is not None:
        with open(attachment, "rb") as f:
            mime_msg.attach(MIMEApplication(
                f.read(),
                Content_Disposition='attachment; filename="%s"' % os.path.basename(attachment),
                Name=os.path.basename(attachment)))

    try:
        smtp = smtplib.SMTP_SSL(host, int(port))
        _, s = smtp.login(user, passwd)
        print(s)
        msg_str = "From: {0}\r\nTo: {1}\r\nSubject: {2}\r\n\r\n{3}".format(from_addr, to_addr, subject, msg)
        d = smtp.sendmail(from_addr, to_addr, mime_msg.as_string())
        if len(d) != 0:
            print(d)
        smtp.quit()
    except (smtplib.SMTPRecipientsRefused, smtplib.SMTPHeloError,
            smtplib.SMTPSenderRefused, smtplib.SMTPDataError) as e:
        print("Server error:", str(e))
    return

if __name__ == "__main__":

    # Send a test message, containing the module as an attachment
    send_mail("mail.conf", "cedar@ironicmtn.com", "njwilson23@gmail.com",
              "testing auto mailer", subject="test", attachment="aspmail.py")

