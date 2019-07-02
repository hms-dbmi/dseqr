#!/bin/bash

read -p "Github username: " user
read -s -p "Github password: " pass

ssh ubuntu@ec2-52-14-124-45.us-east-2.compute.amazonaws.com 'bash -s' < scripts/update-scseq.sh $user $pass