HPCI Shared storage access on Fugaku
************************************
This tutorial explains how to access the HPCI shared storage from Fugaku.
The shared storage system is a large-scale data sharing platform for HPCI users. By
using the shared storage system, HPCI users can quickly and safely share large
amounts of data under one file system across the geographically-dispersed
computational resources of the HPCI.
The HPCI shared storage is different from the default shared storage available
on Fugaku. 
In order to get access to the HPCI shared storage users need to apply for it within 
a project proposal. 


Accessing the shared storage from Fugaku
========================================

#. Issue a client and proxy certificate
   
   Visit the website https://portal.hpci.nii.ac.jp to issue a client and proxy certificate.
   Select Riken CCSE from the dropdown menu. Login with your HPCI user id and password.
   You should see a screen similar to the one shown below.
   Click **Issuing a client certificate**. Enter a passphrase. After issuing the client
   certificate click **Back to Menu**.


   .. image:: img/certificate.png
      :width: 100%

   In the next step click **Storing/Downloading Proxy Certificate**. You should see
   the following screen. Enter the passphrase for the client certificate you entered
   in the previous step. Select the number of hours the certificate should be vaild
   in the **valid for (in hours)** dropdown menu, **168** hours are recommended. 
   Enter a new passphrase for the proxy certificate and click the **Go** button.

   .. image:: img/cert2.png
      :width: 100%

   .. warning::

    The proxy certificate is only valid for a week. You need to repeat the step above
    if the proxy certificate experires and you want to access the shared storage again.

#. Login to the Fugaku Cloud Storage Gateway Node
   
   The Fugaku Cloud Storage Gateway Node is a dedicated node for accessing the
   HPCI shared storage. Only this node can be used to access the shared storage
   from Fugaku. The Cloud Storage Gateway Node has the same compiler and software
   environment as the default Fugaku nodes. Login with the following command:

   .. code-block:: bash

    ssh YourUserId@csgw.fugaku.r-ccs.riken.jp

   After logging in to the Cloud Storage Gateway Node, check the expiration date of the
   proxy certificate with the grid-proxy-info command.
   The following is an example of a case where the proxy certificate has expired. If it has
   expired, you need to obtain a new proxy certificate.

   .. code-block:: bash

    [YourUsername@csgw1 ~]$ grid-proxy-info
       ERROR: Couldnt find a valid proxy.
       globus_sysconfig: Could not find a valid proxy certificate file location
       globus_sysconfig: Error with key filename
       globus_sysconfig: File does not exist: /tmp/x509up_pXXXXX is not a valid file
       Use -debug for further information.

   Use the myproxy-logon command to obtain a proxy certificate as follows, pass your **HPCI id**
   as a parameter to the command:

   .. code-block:: bash
    
    [YourUsername@csgw1 ~]$ myproxy-logon -s portal.hpci.nii.ac.jp -l YourHPCIId -t 168
       Enter MyProxy pass phrase: ******
       A credential has been received for user YourHPCIId in
       /tmp/x509up_XXXXXX.fileXXXXXXX.

   After acquiring the proxy certificate, run the **grid-proxy-info** command again to check
   the validity period.  The validity period will be displayed in the timeleft field

   .. code-block:: bash

    [YourUsername@csgw1 ~]$ grid-proxy-info
     subject : /C=JP/O=NII/OU=HPCI/CN=Hoge%40Foo[hpci000000]/CN=XXXXXXXXXX/CN=XXXXXXXXX/CN=XXXXXX
     issuer : /C=JP/O=NII/OU=HPCI/CN=Hoge%40Foo[hpci000000]/CN=XXXXXXXXXX/CN=XXXXXXXXX/CN=XXXXXX
     identity : /C=JP/O=NII/OU=HPCI/CN=Hoge%40Foo[hpci000000]
     type : RFC 3820 compliant impersonation proxy
     strength : 2048 bits
     path : /tmp/x509up_pXXXX.fileXYZABCD
     timeleft : 23:59:40

   You can mount the shared storage by using the **mount.hpci** command:


   .. code-block:: bash
 
    [u00000@csgw1 ~]$ mount.hpci
      Update proxy certificate for gfarm2fs
      timeleft : 23:53:05
      Mount GfarmFS on /gfarm/hp000000/u000000
      Mount GfarmFS on /gfarm/hp000001/u000000

   .. note::
   
     The mount destination of the shared storage is displayed in the next field of "Mount GfarmFS on".
     Normally, it will be mounted on /gfarm, but if the directory does not exist, it will be mounted on /tmp, as
     in /tmp/hp000000/u000000.

#. Copy files between Fugaku and the HPCI shared storage server
   
   You can use the normal **cp** command to copy files in addition the shared torage provides the **gfpcopy** command to copy 
   multiple files in parallel.
   In the gfpcopy command, the parallelism of the copy is specified by the -j option. The default parallelism
   is 4.

   .. code-block:: bash
 
    gfpcopy -j8 -v my_data /tmp/hp000000/u000000

   .. warning::

    The **gfpcopy** command sometimes hangs and the copy process does not continue.
    The **cp** command is recommended.


   .. note::

    Remove the -v option from the command in order to avoid the output.


   The following table shows other useful commands:

   +--------------------------+------------------------------------------------------------------------------------+
   | Command                  | Explanation                                                                        |
   +==========================+====================================================================================+
   | gfusage                  | The gfusage command outputs the amount of usage and number of files used by users. |
   +--------------------------+------------------------------------------------------------------------------------+
   | gfquota -g ProjectID -H  | The gfquota command outputs the allocated memory and additional information.       |
   +--------------------------+------------------------------------------------------------------------------------+
