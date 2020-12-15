#!groovy
pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
            sh 'cmake -Bbuild -H. && cd build && make'
            }
        }
        stage('Test') {
            steps {
                sh 'cd build && make test'
            }
        }
    }
}
