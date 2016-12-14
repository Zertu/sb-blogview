const express = require('express'),
  router = express.Router(),
  user = require('./user')
  /* GET home page. */
router.get('/', function (req, res, next) {
  res.render('index', {
    title: 'sb-blog'
  });
});
router.get('/article', (req, res) => {

  console.log(req)
  res.render('article', {
    title: 'sb-blog',
    session: req.session.user
  })
})


router.post('/login', (req, res) => {
  let row = user.checklogin(req, res, function (rows) {
    res.render('article', {
      title: 'sb-blog',
      session: req.session.user
    })
  })
})

module.exports = router;