#!/bin/bash
/usr/lib/gvfs/gvfs-fuse-daemon ~/mnt
gvfs-mount SFTP://adminrig@192.168.30.11
gvfs-mount SFTP://sbsuser@192.168.30.12/illumina/runs user=sbsuser password=sbs123


