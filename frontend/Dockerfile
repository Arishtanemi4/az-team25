# Build context is the repo root (see docker-compose.yml), so this only ever touches frontend/.
#
# VITE_*_URL (ragServiceClient.ts, scoringServiceClient.ts, geneServiceClient.ts,
# dataServiceClient.ts, research/client.ts) all default to http://localhost:8000 when unset --
# read at build time and baked into the static bundle, then called from the *browser*, not from
# inside this container. That default is exactly right for docker-compose.yml, which publishes
# backend on host port 8000, so no build args are needed here.
FROM node:20-slim AS build
WORKDIR /app
COPY frontend/package.json frontend/package-lock.json ./
RUN npm ci
COPY frontend .
RUN npm run build

FROM nginx:alpine
COPY --from=build /app/dist /usr/share/nginx/html
COPY frontend/nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 80
