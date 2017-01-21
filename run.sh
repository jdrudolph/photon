#!/bin/bash

/redis-stable/src/redis-server /photon/redis.conf && \
cd /photon/phos && \
celery worker -A app.celery --detach && \
cd / && \
/photon/phos/app.py $@
