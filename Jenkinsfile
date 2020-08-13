def bmarkFile = 'benchmarks.jl'
pipeline {
  agent any
  environment {
    REPO_EXISTS = fileExists "$repo"
  }
  options {
    skipDefaultCheckout true
  }
  triggers {
    GenericTrigger(
     genericVariables: [
        [
            key: 'action', 
            value: '$.action',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '[^(created)]', //Optional, defaults to empty string
            defaultValue: '' //Optional, defaults to empty string
        ],
        [
            key: 'comment',
            value: '$.comment.body',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '', //Optional, defaults to empty string
            defaultValue: '' //Optional, defaults to empty string
        ],
        [
            key: 'org',
            value: '$.organization.login',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '', //Optional, defaults to empty string
            defaultValue: 'ProofOfConceptForJuliSmoothOptimizers' //Optional, defaults to empty string
        ],
        [
            key: 'pullrequest',
            value: '$.issue.number',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '[^0-9]', //Optional, defaults to empty string
            defaultValue: '' //Optional, defaults to empty string
        ],
        [
            key: 'repo',
            value: '$.repository.name',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '', //Optional, defaults to empty string
            defaultValue: '' //Optional, defaults to empty string
        ]
     ],

     causeString: 'Triggered on $comment',

     token: "LDLFactorizations",

     printContributedVariables: true,
     printPostContent: true,

     silentResponse: false,

     regexpFilterText: '$comment',
     regexpFilterExpression: '@JSOBot runbenchmarks'
    )
  }
  stages {
    stage('clone repo') {
      when {
        expression { REPO_EXISTS == 'false' }
      }
      steps {
        sh 'git clone https://${GITHUB_AUTH}@github.com/$org/$repo.git'
      }
    }
    stage('checkout on new branch') {
      steps {
        dir(WORKSPACE + "/$repo") {
          sh '''
          git checkout $BRANCH_NAME
          git clean -fd
          git pull
          git fetch --no-tags origin '+refs/heads/master:refs/remotes/origin/master'
          git checkout -b benchmark
          '''    
        }   
      }
    }
    stage('run benchmarks') {
      steps {
        script {
          def data = env.comment.tokenize(' ')
          if (data.size() > 2) {
            bmarkFile = data.get(2);
          }
        }
        dir(WORKSPACE + "/$repo") {
          sh "set -x"
          sh "qsub -N $repo_$pullrequest -V -cwd -o $HOME/benchmarks/bmark_output.log -e $HOME/benchmarks/bmark_error.log push_benchmarks.sh $bmarkFile"
        }   
      }
    }
  }
  post {
    success {
      echo "SUCCESS!"  
    }
    cleanup {
      dir(WORKSPACE + "/$repo") {
        sh 'printenv'
        sh 'git checkout ' + BRANCH_NAME
        sh '''
        git branch -D benchmark
        git clean -fd
        '''
      }
    }
  }
}
